#include "hierarchical_tree.h"
#include "gestalt_processor2.h"

HierarchicalTree::HierarchicalTree(ScatterPointDataset* data) 
	: TreeCommon(data),
	  fitness_threshold_(0.75),
	  max_level_(4){

}

HierarchicalTree::~HierarchicalTree() {

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
	processor->SetPropertyOn(GestaltProcessor2::PROXIMITY);

	while (current_level < max_level_) {
		float current_max_radius = average_edge_length_ * pow(2, current_level);

		processor->SetDisThreshold(current_max_radius);
		processor->SetThreshold(GestaltProcessor2::SIMILARITY, similarity_threshold_);
		processor->SetThreshold(GestaltProcessor2::PROXIMITY, proximity_threshold_);
		processor->SetData(leaf_nodes_, node_cluster_index_);
		processor->GenerateCluster();

		std::vector< int > new_cluster_index;
		processor->GetCluster(fitness_threshold_, new_cluster_index);

		if (new_cluster_index.size() != 0) {
			for (int i = 0; i < new_cluster_index.size(); ++i)
				for (int j = 0; j < node_cluster_index_.size(); ++j) {
					if (node_cluster_index_[j] == new_cluster_index[i]) node_cluster_index_[j] = new_cluster_index[0];
				}
			for (int i = 1; i < new_cluster_index.size(); ++i)
				for (int j = 0; j < node_cluster_index_.size(); ++j)
					if (node_cluster_index_[j] > new_cluster_index[i]) node_cluster_index_[j] -= 1;
			current_cluster_count -= (new_cluster_index.size() - 1);
		} else {
			current_level++;
		}
	}
}

void HierarchicalTree::GenerateFitValues(int origin_count, int current_count, std::vector< int >& origin_index, std::vector< int >& current_index) {

}