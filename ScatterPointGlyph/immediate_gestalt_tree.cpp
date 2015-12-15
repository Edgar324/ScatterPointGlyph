#include "immediate_gestalt_tree.h"
#include "gestalt_processor2.h"
#include "scatter_point_dataset.h"

ImmediateGestaltTree::ImmediateGestaltTree(ScatterPointDataset* dataset)
	: ImmediateTree(dataset) {
}

ImmediateGestaltTree::~ImmediateGestaltTree() {
}

void ImmediateGestaltTree::run() {
	if (!is_threshold_updated_) return;
	if (root_->linked_nodes.size() != 0) {
		for (int i = 0; i < root_->linked_nodes.size(); ++i) delete root_->linked_nodes[i];
		root_->linked_nodes.clear();
	}
	this->ConstructDirectly();
	this->GenerateCluster();
	is_threshold_updated_ = false;
}

void ImmediateGestaltTree::GenerateCluster() {
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

	while (true) {
		proximity_threshold_ = 1.0 / 3 * max_radius_threshold_;
		processor->SetDisThreshold(max_radius_threshold_);
		processor->SetThreshold(GestaltProcessor2::SIMILARITY, similarity_threshold_);
		processor->SetThreshold(GestaltProcessor2::PROXIMITY, proximity_threshold_);
		processor->SetData(leaf_nodes_, cluster_nodes, node_cluster_index_, node_connecting_status_, dataset_->weights);
		processor->GenerateCluster();

		std::vector< int > new_cluster_index;
		processor->GetCluster(0.0, new_cluster_index);

		if (new_cluster_index.size() != 0) {
			for (int i = 0; i < new_cluster_index.size(); ++i)
				for (int j = 0; j < node_cluster_index_.size(); ++j) {
					if (node_cluster_index_[j] == new_cluster_index[i]) node_cluster_index_[j] = new_cluster_index[0];
				}
			for (int i = new_cluster_index.size() - 1; i >= 1; --i)
				for (int j = 0; j < node_cluster_index_.size(); ++j)
					if (node_cluster_index_[j] > new_cluster_index[i]) node_cluster_index_[j] -= 1;
			current_cluster_count -= (new_cluster_index.size() - 1);

			CBranch* branch_node = new CBranch;
			branch_node->set_level(current_level);
			for (int i = 0; i < new_cluster_index.size(); ++i)
				branch_node->linked_nodes.push_back(cluster_nodes[new_cluster_index[i]]);
			cluster_nodes[new_cluster_index[0]] = branch_node;
			for (int i = new_cluster_index.size() - 1; i >= 1; --i)
				for (int j = new_cluster_index[i]; j < cluster_nodes.size() - 1; ++j)
					cluster_nodes[j] = cluster_nodes[j + 1];
			cluster_nodes.resize(cluster_nodes.size() - new_cluster_index.size() + 1);

			std::vector< float > average_value;
			std::vector< float > average_pos;
			int node_count = 0;
			average_pos.resize(2, 0);
			average_value.resize(dataset_->weights.size(), 0);
			for (int j = 0; j < node_cluster_index_.size(); ++j)
				if (node_cluster_index_[j] == new_cluster_index[0]) {
					for (int k = 0; k < dataset_->weights.size(); ++k)
						average_value[k] += leaf_nodes_[j]->average_values[k];
					average_pos[0] += leaf_nodes_[j]->center_pos[0];
					average_pos[1] += leaf_nodes_[j]->center_pos[1];
					node_count += 1;
				}
			for (int k = 0; k < dataset_->weights.size(); ++k)
				average_value[k] /= node_count;
			average_pos[0] /= node_count;
			average_pos[1] /= node_count;
			branch_node->average_values = average_value;
			branch_node->center_pos = average_pos;
		}
		else {
			break;
		}
	}

	root_->linked_nodes = cluster_nodes;
	root_->set_level(2);
}