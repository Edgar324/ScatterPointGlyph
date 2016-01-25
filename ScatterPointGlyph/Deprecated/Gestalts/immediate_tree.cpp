#include "immediate_tree.h"
#include "scatter_point_dataset.h"
#include "gestalt_processor2.h"

ImmediateTree::ImmediateTree(ScatterPointDataset* data)
	: TreeCommon(data),
	similarity_threshold_(0.3),
	proximity_threshold_(0.1),
	max_radius_threshold_(0.0),
	sample_num_(100),
	is_threshold_updated_(true) {

}

ImmediateTree::~ImmediateTree() {

}

void ImmediateTree::SetSimilarityThreshold(float thre) {
	this->similarity_threshold_ = thre;
	is_threshold_updated_ = true;
}

void ImmediateTree::SetProximityThreshold(float thre) {
	this->proximity_threshold_ = thre;
	is_threshold_updated_ = true;
}

void ImmediateTree::SetRadiusThreshold(float max_radius) {
	this->max_radius_threshold_ = max_radius;
	is_threshold_updated_ = true;
}

void ImmediateTree::SetSampleSize(int num) {
	sample_num_ = num;
	is_threshold_updated_ = true;
}

void ImmediateTree::GetClusterResult(int level, int& cluster_num, std::vector< int >& cluster_index) {
	cluster_index.resize(dataset_->original_point_pos.size());
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i) cluster_index[i] = -1;

	std::vector< CNode* > level_node;
	this->Traverse(root_->level() - 1, level_node);

	cluster_num = level_node.size();

	for (int i = 0; i < level_node.size(); ++i) {
		std::vector< int > point_vec;
		this->Traverse(level_node[i], point_vec);
		for (int j = 0; j < point_vec.size(); ++j) cluster_index[point_vec[j]] = i;
	}
}

void ImmediateTree::run() {
	if (!is_threshold_updated_) return;
	if (root_->linked_nodes.size() != 0) {
		for (int i = 0; i < root_->linked_nodes.size(); ++i) delete root_->linked_nodes[i];
		root_->linked_nodes.clear();
	}
	this->ConstructDirectly();
	this->GenerateCluster();
	is_threshold_updated_ = false;
}

void ImmediateTree::GenerateCluster() {
	/*if (leaf_nodes_.size() == 0) {
		std::cout << "Need data initialization first." << std::endl;
		return;
	}

	node_cluster_index_.resize(leaf_nodes_.size());
	for (int i = 0; i < leaf_nodes_.size(); ++i) node_cluster_index_[i] = i;

	int current_level = 1;
	int current_cluster_count = leaf_nodes_.size();

	GestaltProcessor2* processor = new GestaltProcessor2;
	processor->SetPropertyOn(GestaltProcessor2::SIMILARITY);
	processor->SetPropertyOn(GestaltProcessor2::PROXIMITY);

	std::vector< CNode* > cluster_nodes;
	cluster_nodes.resize(leaf_nodes_.size());
	for (int i = 0; i < leaf_nodes_.size(); ++i) cluster_nodes[i] = leaf_nodes_[i];

	proximity_threshold_ = 1.0 / 3 * max_radius_threshold_;
	processor->SetDisThreshold(max_radius_threshold_);
	processor->SetThreshold(GestaltProcessor2::SIMILARITY, similarity_threshold_);
	processor->SetThreshold(GestaltProcessor2::PROXIMITY, proximity_threshold_);
	processor->SetData(leaf_nodes_, cluster_nodes, node_cluster_index_, node_connecting_status_, dataset_->weights);
	processor->GenerateCluster();
	processor->GetResultLabel(node_cluster_index_);

	std::vector< CBranch* > label_nodes;
	label_nodes.resize(leaf_nodes_.size(), NULL);
	for (int i = 0; i < node_cluster_index_.size(); ++i) {
		int label = node_cluster_index_[i];
		if (label_nodes[label] == NULL) {
			label_nodes[label] = new CBranch;
		}
		label_nodes[label]->linked_nodes.push_back(leaf_nodes_[i]);
	}
	for (int i = 0; i < node_cluster_index_.size(); ++i)
		if (label_nodes[i] != NULL) root_->linked_nodes.push_back(label_nodes[i]);
	root_->set_level(1);*/
}