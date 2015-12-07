#include "hierarchical_tree.h"
#include "gestalt_processor2.h"

HierarchicalTree::HierarchicalTree(ScatterPointDataset* data) 
	: TreeCommon(data),
	  fitness_threshold_(0.75),
	  max_level_(4){
	this->ConstructDirectly();
}

HierarchicalTree::~HierarchicalTree() {

}

void HierarchicalTree::SetVariableWeights(std::vector< float >& weights) {
	var_weights_ = weights;
}

void HierarchicalTree::GenerateCluster(float dis_per_pixel, int min_pixel_radius) {
	if (leaf_nodes_.size() == 0) {
		std::cout << "Need data initialization first." << std::endl;
		return;
	}

	node_cluster_index_.resize(leaf_nodes_.size());
	for (int i = 0; i < leaf_nodes_.size(); ++i) node_cluster_index_[i] = i;

	int current_level = 1;
	int current_cluster_count = leaf_nodes_.size();

	GestaltProcessor2* processor = new GestaltProcessor2;
	processor->SetPropertyOn(GestaltProcessor2::SIMILARITY);
	processor->SetPropertyOff(GestaltProcessor2::PROXIMITY);

	std::vector< CNode* > cluster_nodes;
	cluster_nodes.resize(leaf_nodes_.size());
	for (int i = 0; i < leaf_nodes_.size(); ++i) cluster_nodes[i] = leaf_nodes_[i];
	std::vector< std::vector< bool > > cluster_connecting_status = node_connecting_status_;

	while (current_level < max_level_) {
		float current_max_radius = average_edge_length_ * pow(2, current_level - 1);

		processor->SetDisThreshold(current_max_radius);
		//processor->SetThreshold(GestaltProcessor2::SIMILARITY, similarity_threshold_);
		//processor->SetThreshold(GestaltProcessor2::PROXIMITY, proximity_threshold_);
		processor->SetData(leaf_nodes_, cluster_nodes, node_cluster_index_, node_connecting_status_, var_weights_);
		processor->GenerateCluster();

		std::vector< int > new_cluster_index;
		processor->GetCluster(fitness_threshold_, new_cluster_index);

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
			average_value.resize(var_weights_.size(), 0);
			for (int j = 0; j < node_cluster_index_.size(); ++j) 
				if (node_cluster_index_[j] == new_cluster_index[0]) {
					for (int k = 0; k < var_weights_.size(); ++k)
						average_value[k] += leaf_nodes_[j]->average_values[k];
					average_pos[0] += leaf_nodes_[j]->center_pos[0];
					average_pos[1] += leaf_nodes_[j]->center_pos[1];
					node_count += 1;
				}
			for (int k = 0; k < var_weights_.size(); ++k)
				average_value[k] /= node_count;
			average_pos[0] /= node_count;
			average_pos[1] /= node_count;
			branch_node->average_values = average_value;
			branch_node->center_pos = average_pos;
		} else {
			current_level++;
		}
	}
}