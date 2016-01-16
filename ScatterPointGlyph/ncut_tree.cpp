#include "ncut_tree.h"
#include "scatter_point_dataset.h"
#include <queue>

NCutTree::NCutTree(ScatterPointDataset* data) 
	: TreeCommon(data), un_threshold_(0.2)
{
}

NCutTree::~NCutTree()
{

}

void NCutTree::SetUncertaintyThreshold(float un_threshold)
{
	un_threshold_ = un_threshold;
}

void NCutTree::GetClusterResult(float radius, std::vector< CNode* >& level_nodes)
{

}

void NCutTree::GetClusterResult(float dis_per_piexl, int& cluster_num, std::vector< int >& cluster_index)
{
}

void NCutTree::GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index)
{
}

void NCutTree::run()
{
	if (root_->linked_nodes.size() != 0) {
		for (int i = 0; i < root_->linked_nodes.size(); ++i) delete root_->linked_nodes[i];
		root_->linked_nodes.clear();
	}

	this->ConstructDirectly();

	root_->set_level(0);
	root_->is_expanded = true;
	for (int i = 0; i < root_->linked_nodes.size(); ++i) {
		root_->linked_nodes[i]->set_level(1);
		root_->linked_nodes[i]->parent = root_;
	}
	ProgressNode(root_);


	std::queue< CNode* > node_queue;
	node_queue.push(root_);
	while (node_queue.size() > 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
			bool is_un_fit = false;
			for (int i = 0; i < dataset_->var_num; ++i)
				if (branch == root_ || branch->value_variance[i] > un_threshold_) {
					is_un_fit = true;
					break;
				}
			if (is_un_fit) GenerateCluster(branch);

			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				node_queue.push(branch->linked_nodes[i]);
		}
	}

	AssignColor(root_, 0, 1.0);

	this->InitializeSortingIndex();
}

void NCutTree::GenerateCluster(CBranch* node) 
{
	
}