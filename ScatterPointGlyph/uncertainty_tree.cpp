#include "uncertainty_tree.h"
#include "scatter_point_dataset.h"
#include "multi_label_processor.h"

UncertaintyTree::UncertaintyTree(ScatterPointDataset* data)
	: TreeCommon(data),
	max_radius_threshold_(0.0),
	sample_num_(100),
	is_threshold_updated_(true) {

	processor_ = new MultiLabelProcessor;
}

UncertaintyTree::~UncertaintyTree() {

}

void UncertaintyTree::SetRadiusThreshold(float max_radius) {
	this->max_radius_threshold_ = max_radius;
	is_threshold_updated_ = true;
}

void UncertaintyTree::SetSampleSize(int num) {
	sample_num_ = num;
	is_threshold_updated_ = true;
}

void UncertaintyTree::GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index) {
	cluster_index.clear();

	std::vector< CNode* > level_node;
	this->Traverse(root_->level() - 1, level_node);

	cluster_index.resize(level_node.size());
	for (int i = 0; i < level_node.size(); ++i)
		this->Traverse(level_node[i], cluster_index[i]);
}

void UncertaintyTree::GetClusterResult(float dis_per_pixel, int& cluster_num, std::vector< int >& cluster_index) {
	cluster_index.resize(dataset_->original_point_pos.size());
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i) cluster_index[i] = -1;

	std::vector< CNode* > level_node;
	this->Traverse(root_->level() - 1, level_node);

	cluster_num = level_node.size();

	for (int i = 0; i < level_node.size(); ++i) {
		std::vector< int > point_vec;
		this->Traverse(level_node[i], point_vec);
		for (int j = 0; j < point_vec.size(); ++j) cluster_index[dataset_->sample_index[point_vec[j]]] = i;
	}

	std::cout << "Cluster Count: " << cluster_num << std::endl;
}

std::vector< float >& UncertaintyTree::GetUncertainty() {
	un_value_.resize(dataset_->original_point_pos.size(), 0);

	std::vector< CNode* > level_node;
	this->Traverse(root_->level() - 1, level_node);

	for (int i = 0; i < level_node.size(); ++i) {
		std::vector< int > point_vec;
		this->Traverse(level_node[i], point_vec);
		for (int j = 0; j < point_vec.size(); ++j) un_value_[dataset_->sample_index[point_vec[j]]] = leaf_node_un_[point_vec[j]];
	}

	return un_value_;
}

void UncertaintyTree::run() {
	if (!is_threshold_updated_) return;
	if (root_->linked_nodes.size() != 0) {
		for (int i = 0; i < root_->linked_nodes.size(); ++i) delete root_->linked_nodes[i];
		root_->linked_nodes.clear();
	}

	this->ConstructDirectly();
	for (int i = 0; i < leaf_nodes_.size() - 1; ++i)
		for (int j = i + 1; j < leaf_nodes_.size(); ++j) 
			if (node_connecting_status_[i][j]) {
				float dis = sqrt(pow(leaf_nodes_[i]->center_pos[0] - leaf_nodes_[j]->center_pos[0], 2) + pow(leaf_nodes_[i]->center_pos[1] - leaf_nodes_[j]->center_pos[1], 2));
				if (dis > 0.5 * max_radius_threshold_) {
					node_connecting_status_[i][j] = false;
					node_connecting_status_[j][i] = false;
				}
			}

	edge_weights_.resize(leaf_nodes_.size());
	for (int i = 0; i < leaf_nodes_.size(); ++i)
		edge_weights_[i].resize(leaf_nodes_.size(), 0.5);
	leaf_node_un_.resize(leaf_nodes_.size(), 0.5);

	if (dataset_->weights.size() > 1 ) this->GenerateUncertainty();
	this->GenerateCluster();
	is_threshold_updated_ = false;
}

void UncertaintyTree::GenerateUncertainty() {
	std::vector< std::vector< int > > var_labels;
	var_labels.resize(dataset_->weights.size());

	for (int v = 0; v < dataset_->weights.size(); ++v) {
		std::vector< std::vector< float > > pos;
		std::vector< std::vector< float > > value;
		value.resize(leaf_nodes_.size());
		for (int i = 0; i < leaf_nodes_.size(); ++i) {
			pos.push_back(leaf_nodes_[i]->center_pos);
			value[i].push_back(leaf_nodes_[i]->average_values[v]);
		}

		std::vector< float > weights;
		weights.push_back(1.0);

		processor_->SetLabelEstimationRadius(max_radius_threshold_);
		processor_->SetData(pos, value, weights, node_connecting_status_);
		//processor_->SetSampleNumber(sample_num_);
		processor_->GenerateCluster();
		var_labels[v] = processor_->GetResultLabel();
	}

	for (int i = 0; i < leaf_nodes_.size() - 1; ++i)
		for (int j = i + 1; j < leaf_nodes_.size(); ++j) 
			if (node_connecting_status_[i][j]) {
				float average = 0;
				for (int k = 0; k < dataset_->weights.size(); ++k)
					average += (var_labels[k][i] == var_labels[k][j] ? 0 : 1);
				average /= dataset_->weights.size();
				edge_weights_[i][j] = average;
				edge_weights_[j][i] = edge_weights_[i][j];
			}

	leaf_node_un_.assign(leaf_node_un_.size(), 0);

	std::vector< int > point_edge_num;
	point_edge_num.resize(leaf_nodes_.size(), 0);
	for (int i = 0; i < leaf_nodes_.size() - 1; ++i)
		for (int j = i + 1; j < leaf_nodes_.size(); ++j)
			if (node_connecting_status_[i][j]) {
				leaf_node_un_[i] += edge_weights_[i][j];
				leaf_node_un_[j] += edge_weights_[i][j];
				point_edge_num[i]++;
				point_edge_num[j]++;
			}
	for (int i = 0; i < leaf_nodes_.size(); ++i) 
		if (point_edge_num[i] != 0) leaf_node_un_[i] /= point_edge_num[i];
}

void UncertaintyTree::GenerateCluster() {
	if (leaf_nodes_.size() == 0) {
		std::cout << "Need data initialization first." << std::endl;
		return;
	}

	std::vector< std::vector< float > > pos;
	std::vector< std::vector< float > > value;
	for (int i = 0; i < leaf_nodes_.size(); ++i) {
		pos.push_back(leaf_nodes_[i]->center_pos);
		value.push_back(leaf_nodes_[i]->average_values);
	}

	processor_->SetLabelEstimationRadius(max_radius_threshold_);
	processor_->SetData(pos, value, dataset_->weights, node_connecting_status_);
	processor_->SetEdgeWeights(edge_weights_);
	//processor_->SetSampleNumber(sample_num_);
	processor_->GenerateCluster();
	std::vector< int > node_cluster_index = processor_->GetResultLabel();

	std::vector< CBranch* > label_nodes;
	label_nodes.resize(leaf_nodes_.size(), NULL);
	for (int i = 0; i < node_cluster_index.size(); ++i) {
		int label = node_cluster_index[i];
		if (label < 0) continue;
		if (label_nodes[label] == NULL) {
			label_nodes[label] = new CBranch;
		}
		label_nodes[label]->linked_nodes.push_back(leaf_nodes_[i]);
	}
	for (int i = 0; i < node_cluster_index.size(); ++i)
		if (label_nodes[i] != NULL) root_->linked_nodes.push_back(label_nodes[i]);
	root_->set_level(1);
}