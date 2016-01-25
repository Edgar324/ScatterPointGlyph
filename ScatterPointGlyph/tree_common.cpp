#include "tree_common.h"
#include <queue>
#include <time.h>

#include "scatter_point_dataset.h"
#include "utility.h"

TreeCommon::TreeCommon(ScatterPointDataset* data)
	: dataset_(data),
	min_edge_length_(0),
	max_level_(0) {

	root_ = new CBranch;
	root_->is_expanded = true;
	root_->is_highlighted = false;
}

TreeCommon::~TreeCommon() {

}

void TreeCommon::run() {

}

void TreeCommon::ConstructOnOctree(float thre) {

}

void TreeCommon::ConstructOnKmeans(int basic_cnum) {

}

void TreeCommon::ConstructOnRandomSample(int sample_num) {
	srand((unsigned int)time(0));
	std::vector< bool > is_used;
	is_used.resize(dataset_->point_pos.size(), false);

	root_->linked_nodes.clear();

	float temp_min = 1e10;
	for (int i = 0; i < dataset_->point_pos.size(); ++i)
		if (!is_used[i]) {
			is_used[i] = true;
			CLeaf* temp_leaf = new CLeaf();
			std::vector< int > neighbour_list;
			neighbour_list.push_back(i);
			for (int j = i + 1; j < dataset_->point_pos.size(); ++j){
				float dis = sqrt(pow(dataset_->point_pos[j][0] - dataset_->point_pos[i][0], 2) + pow(dataset_->point_pos[j][1] - dataset_->point_pos[i][1], 2));
				if (dis < 1e-5) {
					neighbour_list.push_back(j);
					is_used[j] = true;
				}
				if (dis < temp_min) temp_min = dis;
			}
			temp_leaf->linked_points = neighbour_list;

			temp_leaf->center_pos.resize(2);
			temp_leaf->center_pos[0] = 0;
			temp_leaf->center_pos[1] = 0;
			temp_leaf->average_values.resize(dataset_->var_weights.size());
			for (int j = 0; j < dataset_->var_weights.size(); ++j)
				temp_leaf->average_values[j] = 0;
			for (int j = 0; j < neighbour_list.size(); ++j) {
				temp_leaf->center_pos[0] += dataset_->point_pos[neighbour_list[j]][0];
				temp_leaf->center_pos[1] += dataset_->point_pos[neighbour_list[j]][1];
				for (int k = 0; k < dataset_->var_weights.size(); ++k)
					temp_leaf->average_values[k] += dataset_->point_values[neighbour_list[j]][k];
			}
			temp_leaf->center_pos[0] /= neighbour_list.size();
			temp_leaf->center_pos[1] /= neighbour_list.size();
			for (int j = 0; j < dataset_->var_weights.size(); ++j)
				temp_leaf->average_values[j] /= neighbour_list.size();

			root_->linked_nodes.push_back(temp_leaf);
			temp_leaf->set_level(0);
		}
}

void TreeCommon::ConstructDirectly() {
	float min_dis_threshold = 1e-5;
	std::vector< bool > is_used;
	is_used.resize(dataset_->point_pos.size(), false);

	root_->linked_nodes.clear();

	float temp_min = 1e10;
	for (int i = 0; i < dataset_->point_pos.size(); ++i) 
		if (!is_used[i]) {
			is_used[i] = true;
			CLeaf* temp_leaf = new CLeaf();
			std::vector< int > neighbour_list;
			neighbour_list.push_back(i);
			for (int j = i + 1; j < dataset_->point_pos.size(); ++j){
				float dis = sqrt(pow(dataset_->point_pos[j][0] - dataset_->point_pos[i][0], 2) + pow(dataset_->point_pos[j][1] - dataset_->point_pos[i][1], 2));
				if (dis < 1e-10) {
					neighbour_list.push_back(j);
					is_used[j] = true;
				}
				if (dis < temp_min) temp_min = dis;
			}
			temp_leaf->linked_points = neighbour_list;

			temp_leaf->center_pos.resize(2);
			temp_leaf->center_pos[0] = 0;
			temp_leaf->center_pos[1] = 0;
			temp_leaf->average_values.resize(dataset_->var_weights.size());
			for (int j = 0; j < dataset_->var_weights.size(); ++j)
				temp_leaf->average_values[j] = 0;
			for (int j = 0; j < neighbour_list.size(); ++j) {
				temp_leaf->center_pos[0] += dataset_->point_pos[neighbour_list[j]][0];
				temp_leaf->center_pos[1] += dataset_->point_pos[neighbour_list[j]][1];
				for (int k = 0; k < dataset_->var_weights.size(); ++k)
					temp_leaf->average_values[k] += dataset_->point_values[neighbour_list[j]][k];
			}
			temp_leaf->center_pos[0] /= neighbour_list.size();
			temp_leaf->center_pos[1] /= neighbour_list.size();
			for (int j = 0; j < dataset_->var_weights.size(); ++j)
				temp_leaf->average_values[j] /= neighbour_list.size();
			temp_leaf->point_count = temp_leaf->linked_points.size();

			root_->linked_nodes.push_back(temp_leaf);
			temp_leaf->set_level(0);
		}

	min_edge_length_ = 1e10;
	for (int i = 0; i < root_->linked_nodes.size() - 1; ++i)
		for (int j = i + 1; j < root_->linked_nodes.size(); ++j) {
			float temp_dis;
			temp_dis = sqrt(pow(root_->linked_nodes[i]->center_pos[0] - root_->linked_nodes[j]->center_pos[0], 2)
				+ pow(root_->linked_nodes[i]->center_pos[1] - root_->linked_nodes[j]->center_pos[1], 2));
			if (temp_dis < min_edge_length_) min_edge_length_ = temp_dis;
		}
	/*if (dataset_->is_structured_data) {
		node_connecting_status_.resize(dataset_->point_pos.size());
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			node_connecting_status_[i].resize(dataset_->point_pos.size());
			node_connecting_status_[i].assign(node_connecting_status_[i].size(), false);
		}
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			int x = dataset_->sample_index[i] % dataset_->w;
			int y = dataset_->sample_index[i] / dataset_->w;
			if (x == dataset_->w - 1 || y == dataset_->h - 1) continue;
			int pre_one = -1;
			if (dataset_->node_sample_map.find(y * dataset_->w + x + 1) != dataset_->node_sample_map.end())
				pre_one = dataset_->node_sample_map[y * dataset_->w + x + 1];
			if (pre_one == -1) continue;
			int pre_two = -1;
			if (dataset_->node_sample_map.find((y + 1) * dataset_->w + x) != dataset_->node_sample_map.end())
				pre_two = dataset_->node_sample_map[(y + 1) * dataset_->w + x];
			if (pre_two == -1) continue;
			node_connecting_status_[pre_one][i] = true;
			node_connecting_status_[i][pre_one] = true;
			node_connecting_status_[pre_two][i] = true;
			node_connecting_status_[i][pre_two] = true;
		}

		this->min_edge_length_ = sqrt(pow(leaf_nodes_[0]->center_pos[0] - leaf_nodes_[1]->center_pos[0], 2)
			+ pow(leaf_nodes_[0]->center_pos[1] - leaf_nodes_[1]->center_pos[1], 2));
	}
	else {
		
	}*/
}

void TreeCommon::GetClusterResult(int level, std::vector< CNode* >& level_nodes) {
	this->Traverse(level, level_nodes);
}

void TreeCommon::GetClusterResult(int level, int& cluster_num, std::vector< int >& cluster_index) {
	cluster_index.resize(dataset_->point_num);
	for (int i = 0; i < dataset_->point_num; ++i) cluster_index[i] = -1;

	std::vector< CNode* > level_node;
	this->Traverse(level, level_node);

	cluster_num = level_node.size();

	for (int i = 0; i < level_node.size(); ++i) {
		std::vector< int > point_vec;
		this->Traverse(level_node[i], point_vec);
		for (int j = 0; j < point_vec.size(); ++j) cluster_index[point_vec[j]] = i;
	}

#ifdef DEBUG_ON
	std::cout << "Cluster Count: " << cluster_num << std::endl;
#endif
}

void TreeCommon::GenerateCluster(CBranch* node /* = NULL */) {

}

void TreeCommon::Traverse(int level, std::vector< CNode* >& nodes) {
	nodes.clear();
	Traverse(level, root_, nodes);
}

void TreeCommon::Traverse(int level, CNode* root, std::vector< CNode* >& nodes) {
	root->is_expanded = true;

	if (root->level() == level) {
		nodes.push_back(root);
		root->is_expanded = false;
	}
	else {
		if (root->type() == CNode::BRANCH && root->level() < level) {
			CBranch* branch = dynamic_cast<CBranch*>(root);
			bool is_all_leaf = true;
			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				if (branch->linked_nodes[i]->type() != CNode::LEAF) {
					is_all_leaf = false;
					break;
				}
			if (is_all_leaf) {
				nodes.push_back(root);
				root->is_expanded = false;
			}
			else {
				branch->is_highlighted = false;
				for (int i = 0; i < branch->linked_nodes.size(); ++i)
					Traverse(level, branch->linked_nodes[branch->sorting_index[i]], nodes);
			}
		}
	}
}

void TreeCommon::Traverse(CNode* node, std::vector< int >& linked_points) {
	std::queue< CNode* > node_queue;
	node_queue.push(node);

	while (node_queue.size() != 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast< CBranch* >(temp_node);
			if (branch != NULL) {
				for (int i = 0; i < branch->linked_nodes.size(); ++i)
					node_queue.push(branch->linked_nodes[i]);
			}
		}
		else if (temp_node->type() == CNode::LEAF) {
			CLeaf* leaf = dynamic_cast< CLeaf* >(temp_node);
			if (leaf != NULL) {
				for (int i = 0; i < leaf->linked_points.size(); ++i)
					linked_points.push_back(leaf->linked_points[i]);
			}
		}
	}
}

void TreeCommon::AssignColor(CNode* node, float hstart, float hend, float factor, bool perm, bool rev) {
	node->color.setHsl((hstart + hend) / 2 * 255, 0.6 * 255, 0.7 * 255);
	node->hstart = hstart;
	node->hend = hend;
	if (node->type() == CNode::LEAF) return;

	CBranch* branch = (CBranch*)node;
	int count = branch->linked_nodes.size();

	float step = (hend - hstart) / count;
	std::vector< float > temp_start, temp_end;
	temp_start.resize(count);
	temp_end.resize(count);

	for (int i = 0; i < count; ++i){
		temp_start[i] = hstart + i * step;
		temp_end[i] = hstart + (i + 1) * step;
	}

	if (perm) {
		for (int i = 0; i < count / 2; ++i) {
			float temp = temp_start[i];
			temp_start[i] = temp_start[count - i - 1];
			temp_start[count - i - 1] = temp;

			temp = temp_end[i];
			temp_end[i] = temp_end[count - i - 1];
			temp_end[count - i - 1] = temp;
		}
	}

	if (rev) {
		int temp_count = (count - 1) / 4;
		for (int i = 0; i < temp_count; ++i) {
			float temp = temp_start[i * 2];
			temp_start[i * 2] = temp_start[temp_count * 4 - i * 2];
			temp_start[temp_count * 4 - i * 2] = temp;

			temp = temp_end[i * 2];
			temp_end[i * 2] = temp_end[temp_count * 4 - i * 2];
			temp_end[temp_count * 4 - i * 2] = temp;
		}
	}

	for (int i = 0; i < branch->linked_nodes.size(); ++i) {
		AssignColor(branch->linked_nodes[i], temp_start[i] + (1.0 - factor) / 2.0 * step, temp_end[i] - (1.0 - factor) / 2.0 * step, factor, perm, rev);
	}
}

void TreeCommon::InitializeSortingIndex() {
	std::queue< CNode* > node_queue;
	node_queue.push(root_);

	while (node_queue.size() != 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
			if (branch != NULL) {
				branch->sorting_index.resize(branch->linked_nodes.size());
				for (int i = 0; i < branch->sorting_index.size(); ++i)
					branch->sorting_index[i] = i;

				for (int i = 0; i < branch->linked_nodes.size(); ++i)
					if (branch->linked_nodes[i]->type() == CNode::BRANCH) node_queue.push(branch->linked_nodes[i]);
			}
		}
	}
}

void TreeCommon::SortTree(std::vector< int >& node_ids) {
	int node_count = 0;
	this->SortNode(root_, node_ids, node_count);
}

int TreeCommon::SortNode(CNode* node, std::vector< int >& node_ids, int& node_count) {
	if (node_count > node_ids.size()) return 9999;

	// find all the nodes in the linked_nodes that exist in the node_index
	if (node->type() == CNode::BRANCH) {
		CBranch* branch = dynamic_cast<CBranch*>(node);
		for (int i = 0; i < node_ids.size(); ++i)
			if (node_ids[i] == node->id()) return i;

		std::vector< bool > is_node_exist;
		std::vector< int > new_sequence;
		std::vector< int > node_min_index;
		node_min_index.resize(branch->linked_nodes.size(), 9999);
		is_node_exist.resize(branch->linked_nodes.size(), false);

		for (int i = 0; i < branch->linked_nodes.size(); ++i)
			if (node_count < node_ids.size()) node_min_index[i] = SortNode(branch->linked_nodes[i], node_ids, node_count);
		for (int i = 0; i < node_min_index.size(); ++i)
			if (node_min_index[i] < 9999) is_node_exist[i] = true;

		new_sequence.resize(branch->linked_nodes.size());
		for (int i = 0; i < branch->linked_nodes.size(); ++i) new_sequence[i] = i;
		Utility::Sort(node_min_index, new_sequence);

		branch->sorting_index.clear();
		for (int i = 0; i < new_sequence.size(); ++i)
			if (node_min_index[i] < 9999) 
				branch->sorting_index.push_back(new_sequence[i]);
			else 
				break;

		for (int i = 0; i < is_node_exist.size(); ++i)
			if (!is_node_exist[i]) branch->sorting_index.push_back(i);

		return node_min_index[0];
	}

	return 9999;
}

void TreeCommon::ProgressNode(CNode* node) {
	if (node == NULL) return;

	// update the statistics
	std::vector< int > point_index;
	this->Traverse(node, point_index);

	std::vector< float > average, variance, center_pos;
	average.resize(dataset_->var_num, 0);
	variance.resize(dataset_->var_num, 0);
	center_pos.resize(2, 0);
	for (int i = 0; i < point_index.size(); ++i) {
		center_pos[0] += dataset_->point_pos[point_index[i]][0];
		center_pos[1] += dataset_->point_pos[point_index[i]][1];

		for (int j = 0; j < dataset_->var_num; ++j)
			average[j] += dataset_->point_values[point_index[i]][j];
	}
	for (int i = 0; i < dataset_->var_num; ++i)
		average[i] /= point_index.size();
	center_pos[0] /= point_index.size();
	center_pos[1] /= point_index.size();

	for (int i = 0; i < point_index.size(); ++i)
		for (int j = 0; j < dataset_->var_num; ++j)
			variance[j] += pow(dataset_->point_values[point_index[i]][j] - average[j], 2);
	for (int i = 0; i < dataset_->var_num; ++i)
		variance[i] = sqrt(variance[i] / point_index.size());
	node->average_values = average;
	node->value_variance = variance;
	node->point_count = point_index.size();
	node->center_pos = center_pos;
}

void TreeCommon::ResetLevel(CNode* node, int level) {
	if (node == NULL) return;
	node->set_level(level);

	if (node->type() == CNode::BRANCH) {
		CBranch* branch = dynamic_cast<CBranch*>(node);
		if (branch != NULL) {
			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				ResetLevel(branch->linked_nodes[i], level + 1);
		}
	}
}