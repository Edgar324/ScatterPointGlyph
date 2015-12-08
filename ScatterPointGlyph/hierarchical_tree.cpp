#include "hierarchical_tree.h"
#include <queue>
#include "gestalt_processor2.h"
#include "scatter_point_dataset.h"

HierarchicalTree::HierarchicalTree(ScatterPointDataset* data) 
	: TreeCommon(data),
	similarity_threshold_(0.3),
	proximity_threshold_(0.1),
	fitness_threshold_(0.5),
	min_pixel_radius_(30),
	max_level_(-1){

	this->ConstructDirectly();
}

HierarchicalTree::~HierarchicalTree() {

}

void HierarchicalTree::SetMaxLevel(int level) {
	this->max_level_ = level;
}

void HierarchicalTree::SetSimilarityThreshold(float thre) {
	this->similarity_threshold_ = thre;
}

void HierarchicalTree::SetProximityThreshold(float thre) {
	this->proximity_threshold_ = thre;
}

void HierarchicalTree::SetFitnessThreshold(float thre) {
	this->fitness_threshold_ = thre;
}

void HierarchicalTree::SetRadiusParameters(float min_pixel_radius) {
	this->min_pixel_radius_ = min_pixel_radius;
}

void HierarchicalTree::GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index) {
	cluster_index.clear();

	int level = (int)(log(dis_per_pixel * min_pixel_radius_ / this->average_edge_length_) / log(2.0) + 1);
	if (level < 0) level = 0;
	if (level > max_level_) level = max_level_;

	std::vector< CNode* > level_node;
	this->Traverse(level, level_node);

	cluster_index.resize(level_node.size());
	for (int i = 0; i < level_node.size(); ++i)
		this->Traverse(level_node[i], cluster_index[i]);
}

void HierarchicalTree::GetClusterResult(float dis_per_pixel, int& cluster_num, std::vector< int >& cluster_index) {
	int level = (int)(log(dis_per_pixel * min_pixel_radius_ / this->average_edge_length_) / log(2.0) + 1);
	if (level < 0) level = 0;
	if (level > max_level_) level = max_level_;

	cluster_index.resize(dataset_->point_pos.size());

	std::vector< CNode* > level_node;
	this->Traverse(level, level_node);

	cluster_num = level_node.size();

	for (int i = 0; i < level_node.size(); ++i) {
		std::vector< int > point_vec;
		this->Traverse(level_node[i], point_vec);
		for (int j = 0; j < point_vec.size(); ++j) cluster_index[point_vec[j]] = i;
	}
}

void HierarchicalTree::run() {
	if (root_->linked_nodes.size() == 0) this->GenerateCluster(min_pixel_radius_);
}

void HierarchicalTree::GenerateCluster(int min_pixel_radius) {
	if (leaf_nodes_.size() == 0) {
		std::cout << "Need data initialization first." << std::endl;
		return;
	}

	if (max_level_ == -1) {
		max_level_ = (int)(log(0.4 / average_edge_length_) / log(2.0)) + 1;
	}

	node_cluster_index_.resize(leaf_nodes_.size());
	for (int i = 0; i < leaf_nodes_.size(); ++i) node_cluster_index_[i] = i;

	int current_level = 1;
	int current_cluster_count = leaf_nodes_.size();

	GestaltProcessor2* processor = new GestaltProcessor2;
	processor->SetPropertyOff(GestaltProcessor2::SIMILARITY);
	processor->SetPropertyOn(GestaltProcessor2::PROXIMITY);

	std::vector< CNode* > cluster_nodes;
	cluster_nodes.resize(leaf_nodes_.size());
	for (int i = 0; i < leaf_nodes_.size(); ++i) cluster_nodes[i] = leaf_nodes_[i];
	std::vector< std::vector< bool > > cluster_connecting_status = node_connecting_status_;

	while (current_level < max_level_) {
		float current_max_radius = this->average_edge_length_ * pow(2, current_level - 1);

		processor->SetDisThreshold(current_max_radius);
		processor->SetThreshold(GestaltProcessor2::SIMILARITY, similarity_threshold_);
		processor->SetThreshold(GestaltProcessor2::PROXIMITY, proximity_threshold_);
		processor->SetData(leaf_nodes_, cluster_nodes, node_cluster_index_, node_connecting_status_, dataset_->weights);
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
		} else {
			current_level++;
		}
	}

	root_->linked_nodes = cluster_nodes;
	root_->set_level(max_level_);
}