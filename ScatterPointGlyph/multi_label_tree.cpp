#include "multi_label_tree.h"
#include "scatter_point_dataset.h"
#include "multi_label_processor.h"
#include <queue>
#include "tour_path_generator.h"
#include "utility.h"

MultiLabelTree::MultiLabelTree(ScatterPointDataset* data)
	: TreeCommon(data),
	max_radius_threshold_(0.5),
	un_threshold_(0.1) {

	processor_ = new MultiLabelProcessor;

	common_parent_node_ = NULL;
}

MultiLabelTree::~MultiLabelTree() {

}

void MultiLabelTree::SetRadiusThreshold(float max_radius) {
	this->max_radius_threshold_ = max_radius;
}

void MultiLabelTree::SetUncertaintyThreshold(float un_threshold) {
	this->un_threshold_ = un_threshold;
}

int MultiLabelTree::GetRadiusLevel(float radius) {
	int level = 0;
	float temp_radius = max_radius_threshold_;
	while (temp_radius / factor_ > radius * 3) {
		temp_radius /= factor_;
		level++;
	}
	if (level > max_level_) level = max_level_;

	return level;
}

void MultiLabelTree::GenerateClusters() {
	id_node_map_.insert(std::map< int, CNode* >::value_type(root_->id(), root_));
	for (int i = 0; i < root_->linked_nodes.size(); ++i) {
		id_node_map_.insert(std::map< int, CNode* >::value_type(root_->linked_nodes[i]->id(), root_->linked_nodes[i]));
	}
	root_->radius = max_radius_threshold_;

	std::queue< CNode* > node_queue;
	node_queue.push(root_);
	while (node_queue.size() > 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
			if (branch == NULL || branch->radius / factor_ < this->min_edge_length_) continue;
			float accu_un = 0.0;
			for (int i = 0; i < dataset_->var_num; ++i) {
				accu_un += branch->variable_variances[i] * dataset_->var_weights[i];
			}

			/*bool is_un_fit = false;
			for (int i = 0; i < dataset_->var_num; ++i)
				if (branch == root_ || branch->value_variance[i] > un_threshold_) {
					is_un_fit = true;
					break;
				}
			if (is_un_fit) SplitNode(branch);*/
			if (branch == root_ || accu_un > un_threshold_) SplitNode(branch);

			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				node_queue.push(branch->linked_nodes[i]);
		}
	}
}

void MultiLabelTree::GenerateSegmentUncertainty(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status, std::vector< std::vector< float > >& edge_weight) {
	std::vector< std::vector< int > > var_labels;
	var_labels.resize(dataset_->var_weights.size());

	for (int v = 0; v < dataset_->var_weights.size(); ++v) {
		std::vector< std::vector< float > > pos;
		std::vector< std::vector< float > > value;
		value.resize(nodes.size());
		for (int i = 0; i < nodes.size(); ++i) {
			pos.push_back(nodes[i]->center_pos);
			value[i].push_back(nodes[i]->average_values[v]);
		}

		std::vector< float > weights;
		weights.push_back(1.0);

		processor_->SetLabelEstimationRadius(max_radius_threshold_);
		processor_->SetData(pos, value, weights, connecting_status);
		processor_->GenerateCluster();
		var_labels[v] = processor_->GetResultLabel();
	}

	edge_weight.resize(nodes.size());
	for (int i = 0; i < nodes.size(); ++i) {
		edge_weight[i].resize(nodes.size());
		edge_weight[i].assign(nodes.size(), 0);
	}

	for (int i = 0; i < nodes.size() - 1; ++i)
		for (int j = i + 1; j < nodes.size(); ++j)
			if (connecting_status[i][j]) {
				float average = 0;
				for (int k = 0; k < dataset_->var_weights.size(); ++k)
					average += (var_labels[k][i] == var_labels[k][j] ? 0 : 1);
				average /= dataset_->var_weights.size();
				edge_weight[i][j] = average;
				edge_weight[j][i] = edge_weight[i][j];
			}
}

void MultiLabelTree::SplitNode(CBranch* node) {
	if (node == NULL || node->linked_nodes.size() == 0) {
		std::cout << "Need data initialization first." << std::endl;
		return;
	}

	std::vector< std::vector< float > > pos;
	std::vector< std::vector< float > > value;
	for (int i = 0; i < node->linked_nodes.size(); ++i) {
		pos.push_back(node->linked_nodes[i]->center_pos);
		value.push_back(node->linked_nodes[i]->average_values);
	}

	std::vector< std::vector< bool > > connecting_status;
	float temp_min_length;
	Utility::VtkTriangulation(node->linked_nodes, connecting_status, temp_min_length);
	for (int i = 0; i < node->linked_nodes.size() - 1; ++i)
		for (int j = i + 1; j < node->linked_nodes.size(); ++j)
			if (connecting_status[i][j]) {
				float dis = sqrt(pow(node->linked_nodes[i]->center_pos[0] - node->linked_nodes[j]->center_pos[0], 2) + pow(node->linked_nodes[i]->center_pos[1] - node->linked_nodes[j]->center_pos[1], 2));
				if (dis > 0.5 * max_radius_threshold_) {
					connecting_status[i][j] = false;
					connecting_status[j][i] = false;
				}
			}

	std::vector< std::vector< float > > edge_weights;
	edge_weights.resize(node->linked_nodes.size());
	for (int i = 0; i < node->linked_nodes.size(); ++i)
		edge_weights[i].resize(node->linked_nodes.size(), 0.5);

	processor_->SetLabelEstimationRadius(node->radius / factor_);
	processor_->SetData(pos, value, dataset_->var_weights, connecting_status);
	//this->GenerateSegmentUncertainty(node->linked_nodes, connecting_status, edge_weights);
	//processor_->SetEdgeWeights(edge_weights);
	processor_->GenerateCluster();
	std::vector< int > node_cluster_index = processor_->GetResultLabel();

	std::vector< CBranch* > label_nodes;
	label_nodes.resize(node->linked_nodes.size(), NULL);
	for (int i = 0; i < node_cluster_index.size(); ++i) {
		int label = node_cluster_index[i];
		if (label < 0) continue;
		if (label_nodes[label] == NULL) {
			label_nodes[label] = new CBranch;
			label_nodes[label]->set_level(node->level() + 1);
			label_nodes[label]->radius = node->radius / factor_;
			label_nodes[label]->parent = node;

			id_node_map_.insert(std::map< int, CNode* >::value_type(label_nodes[label]->id(), label_nodes[label]));
		}
		label_nodes[label]->linked_nodes.push_back(node->linked_nodes[i]);
		node->linked_nodes[i]->set_level(node->level() + 2);
		node->linked_nodes[i]->parent = label_nodes[label];
	}

	common_parent_node_ = node;

	node->linked_nodes.clear();

	std::vector< CNode* > linked_nodes;
	for (int i = 0; i < node_cluster_index.size(); ++i)
		if (label_nodes[i] != NULL) {
			node->linked_nodes.push_back(label_nodes[i]);
			linked_nodes.push_back(label_nodes[i]);
		}
	for (int i = 0; i < linked_nodes.size(); ++i)
		ProgressNodeAndParent(linked_nodes[i]);
	std::vector< int > tour_list;
	TourPathGenerator::GenerateRoundPath(linked_nodes, tour_list);
	node->linked_nodes.clear();
	for (int i = 0; i < tour_list.size(); ++i)
		node->linked_nodes.push_back(linked_nodes[tour_list[i]]);

	/*for (int i = 0; i < node_cluster_index.size(); ++i) 
		if (label_nodes[i] != NULL) {
			node->linked_nodes.push_back(label_nodes[i]);
			ProgressNodeAndParent(label_nodes[i]);
		}*/
}