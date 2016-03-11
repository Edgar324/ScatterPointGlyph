#include "tree_common.h"
#include <queue>
#include <time.h>

#include "scatter_point_dataset.h"
#include "scatter_grid_dataset.h"
#include "utility.h"

TreeCommon::TreeCommon(ScatterPointDataset* data)
	: dataset_(data),
	root_(NULL),
	min_edge_length_(0),
	data_dis_scale_(1.0),
	max_level_(-1),
	tree_mode_(EXPLORATION_MODE) {
}

TreeCommon::~TreeCommon() {

}

void TreeCommon::SetTreeMode(TreeMode mode) {
	tree_mode_ = mode;
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
			temp_leaf->parent = root_;

			root_->linked_nodes.push_back(temp_leaf);
			temp_leaf->set_level(1);
		}

	min_edge_length_ = 1e10;
	for (int i = 0; i < root_->linked_nodes.size() - 1; ++i)
		for (int j = i + 1; j < root_->linked_nodes.size(); ++j) {
			float temp_dis;
			temp_dis = sqrt(pow(root_->linked_nodes[i]->center_pos[0] - root_->linked_nodes[j]->center_pos[0], 2)
				+ pow(root_->linked_nodes[i]->center_pos[1] - root_->linked_nodes[j]->center_pos[1], 2));
			if (temp_dis < min_edge_length_) min_edge_length_ = temp_dis;
		}

	this->ProgressNodeAndParentData(root_);
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

void TreeCommon::run() {
	if (root_ != NULL) delete root_;
	root_ = new CBranch;
	root_->set_level(0);
	root_->is_expanded = true;
	root_->is_highlighted = false;

	this->ConstructDirectly();

    id_node_map_.insert(std::map< int, CNode* >::value_type(root_->id(), root_));
	for (int i = 0; i < root_->linked_nodes.size(); ++i) {
		id_node_map_.insert(std::map< int, CNode* >::value_type(root_->linked_nodes[i]->id(), root_->linked_nodes[i]));
	}

    if (dataset_->type() == ScatterPointDataset::POINT_DATA) {
        Utility::VtkTriangulation(root_->linked_nodes, node_connecting_status_, min_edge_length_);
    }
    else {
        ScatterGridDataset* grid_data = dynamic_cast<ScatterGridDataset*>(dataset_);
        Utility::GridConnection(root_->linked_nodes, grid_data->w, grid_data->h, node_connecting_status_, min_edge_length_);

        grid_id_seq_map_.clear();
        for (int i = 0; i < root_->linked_nodes.size(); ++i)
            grid_id_seq_map_.insert(std::map<int, int>::value_type(root_->linked_nodes[i]->id(), i));
    }

	if (tree_mode_ == EXPLORATION_MODE)
		this->BeginClustering();
	else
		this->GenerateClusters();

	this->ResetSortingIndex(root_);

	this->ResetLevel(root_, 0);

	this->AssignLeafLevel(root_, max_level_);

	this->AssignColor(root_, 0, 1.0);
}

void TreeCommon::Traverse(int level, std::vector< CNode* >& nodes) {
    nodes.clear();
    TraverseLevelNodes(level, root_, nodes);
	
	/*if (tree_mode_ == VIEWING_MODE)
		TraverseLevelNodes(level, root_, nodes);
	else
		TraversAllNodes(root_, nodes);*/
}

void TreeCommon::TraversAllNodes(CNode* root_node, std::vector< CNode* >& nodes) {
	if (root_node != NULL && root_node->type() == CNode::BRANCH) {
		CBranch* branch = dynamic_cast<CBranch*>(root_node);
		branch->is_expanded = true;
		branch->is_highlighted = false;

		bool is_expandable = true;
		for (int i = 0; i < branch->linked_nodes.size(); ++i) 
			if (branch->linked_nodes[i]->type() == CNode::LEAF) {
				is_expandable = false;
				break;
			}
		if (!is_expandable) {
			branch->is_expanded = false;
			nodes.push_back(branch);
		} else {
			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				TraversAllNodes(branch->linked_nodes[branch->sorting_index[i]], nodes);
		}
	}
}

void TreeCommon::TraverseLevelNodes(int level, CNode* root_node, std::vector< CNode* >& nodes) {
	root_node->is_expanded = true;

	if (root_node->level() == level) {
		nodes.push_back(root_node);
		root_node->is_expanded = false;
	}
	else {
		if (root_node->type() == CNode::BRANCH && root_node->level() < level) {
			CBranch* branch = dynamic_cast<CBranch*>(root_node);
			bool is_all_leaf = true;
			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				if (branch->linked_nodes[i]->type() != CNode::LEAF) {
					is_all_leaf = false;
					break;
				}
			if (is_all_leaf && level != max_level_) {
				nodes.push_back(root_node);
				root_node->is_expanded = false;
			}
			else {
				branch->is_highlighted = false;
				for (int i = 0; i < branch->linked_nodes.size(); ++i)
					TraverseLevelNodes(level, branch->linked_nodes[branch->sorting_index[i]], nodes);
			}
		} else if (root_node->type() == CNode::LEAF && root_node->level() < level) {
			nodes.push_back(root_node);
			root_node->is_expanded = false;
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

void TreeCommon::AssignLeafLevel(CNode* node, int level) {
	std::queue< CNode* > node_queue;
	node_queue.push(node);

	while (node_queue.size() != 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
			if (branch != NULL) {
				for (int i = 0; i < branch->linked_nodes.size(); ++i)
					node_queue.push(branch->linked_nodes[i]);
			}
		}
		else if (temp_node->type() == CNode::LEAF) {
			temp_node->set_level(level);
		}
	}
}

void TreeCommon::AssignColor(CNode* node, float hstart, float hend, float factor, bool perm, bool rev) {
	node->color.setHsl((hstart + hend) / 2 * 255, 0.6 * 255, 0.7 * 255);
	node->hstart = hstart;
	node->hend = hend;

	int count = 1;
	if (node->type() != CNode::LEAF) {
		CBranch* branch = (CBranch*)node;
		count = branch->linked_nodes.size();
	}

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

	if (node->type() != CNode::LEAF) {
		CBranch* branch = (CBranch*)node;
		for (int i = 0; i < branch->linked_nodes.size(); ++i) {
			AssignColor(branch->linked_nodes[i], temp_start[i] + (1.0 - factor) / 2.0 * step, temp_end[i] - (1.0 - factor) / 2.0 * step, factor, perm, rev);
		}
	}
}

void TreeCommon::ResetSortingIndex(CNode* node) {
	std::queue< CNode* > node_queue;
	node_queue.push(node);

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

void TreeCommon::ResetLevel(CNode* node, int level) {
	if (node == NULL) return;
	node->set_level(level);
	if (level > max_level_) max_level_ = level;

	if (node->type() == CNode::BRANCH) {
		CBranch* branch = dynamic_cast<CBranch*>(node);
		if (branch != NULL) {
			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				ResetLevel(branch->linked_nodes[i], level + 1);
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

void TreeCommon::SplitCluster(int cluster_index) {
	std::map< int, CNode* >::iterator node_iter = id_node_map_.find(cluster_index);
	if (node_iter == id_node_map_.end()) {
		std::cout << "Cluster not found!" << std::endl;
		return;
	}

	if (node_iter->second->type() != CNode::BRANCH) {
		std::cout << "Leaf or unknown cluster is not allowed to be split!" << std::endl;
		return;
	}

	CBranch* branch = (CBranch*)(node_iter->second);
	bool is_child_leaf = true;
	for (int i = 0; i < branch->linked_nodes.size(); ++i)
		if (branch->linked_nodes[i]->type() != CNode::LEAF) {
			is_child_leaf = false;
			break;
		}
	if (is_child_leaf) this->SplitNode(branch);

	ResetLevel(branch, branch->level());
	ResetSortingIndex(branch);
	AssignColor(branch, branch->hstart, branch->hend);

	/*CBranch* temp_parent = branch->parent;
	if (temp_parent == NULL) return;

	int node_index = -1;
	for (int j = 0; j < temp_parent->linked_nodes.size(); ++j)
		if (temp_parent->linked_nodes[j] == branch) {
			node_index = j;
			break;
		}
	if (node_index != -1) {
		int original_size = temp_parent->linked_nodes.size();
		temp_parent->linked_nodes.resize(temp_parent->linked_nodes.size() + branch->linked_nodes.size() - 1);
		for (int i = original_size; i > node_index; --i)
			temp_parent->linked_nodes[i + branch->linked_nodes.size() - 2] = temp_parent->linked_nodes[i];
		for (int i = 0; i < branch->linked_nodes.size(); ++i) {
			temp_parent->linked_nodes[i + node_index] = branch->linked_nodes[i];
			branch->linked_nodes[i]->parent = temp_parent;
		}

		temp_parent->sorting_index.resize(temp_parent->linked_nodes.size());
		for (int i = 0; i < temp_parent->sorting_index.size(); ++i)
			temp_parent->sorting_index[i] = i;

		ResetLevel(temp_parent, temp_parent->level());
		AssignColor(temp_parent, temp_parent->hstart, temp_parent->hend);
	}*/
}

void TreeCommon::MergeClusters(std::vector< int >& cluster_index) {
	std::vector< CNode* > cluster_nodes;
	for (int i = 0; i < cluster_index.size(); ++i) {
		std::map< int, CNode* >::iterator node_iter = id_node_map_.find(cluster_index[i]);
		if (node_iter == id_node_map_.end()) {
			std::cout << "Cluster not found!" << std::endl;
			return;
		}
		cluster_nodes.push_back(node_iter->second);
	}

	bool is_all_branch = true;
	for (int i = 0; i < cluster_nodes.size(); ++i) {
		if (cluster_nodes[i]->type() != CNode::BRANCH) {
			is_all_branch = false;
			break;
		}
	}

	if (is_all_branch) {
		// find the maximum level node
		int min_level_index = 10000;
		CNode* min_level_node = NULL;
		for (int i = 0; i < cluster_nodes.size(); ++i)
			if (cluster_nodes[i]->level() < min_level_index) {
				min_level_index = cluster_nodes[i]->level();
				min_level_node = cluster_nodes[i];
			}
		if (min_level_node != NULL) {
			CBranch* parent_node = min_level_node->parent;
			CBranch* new_branch = new CBranch;
			new_branch->set_level(min_level_node->level());
			new_branch->radius = min_level_node->radius;
			new_branch->parent = parent_node;
			parent_node->linked_nodes.push_back(new_branch);
			id_node_map_.insert(std::map< int, CNode* >::value_type(new_branch->id(), new_branch));

			this->FindCommonParent(root_, cluster_index);

			for (int i = 0; i < cluster_nodes.size(); ++i) {
				RemoveChildNode(cluster_nodes[i], true);
				CBranch* branch = dynamic_cast<CBranch*>(cluster_nodes[i]);
				for (int j = 0; j < branch->linked_nodes.size(); ++j)
					new_branch->linked_nodes.push_back(branch->linked_nodes[j]);
			}

			ResetSortingIndex(new_branch);
			UpdateChildLevel(new_branch->parent);
			ProgressNodeAndParentData(new_branch);
			AssignColor(common_parent_node_, common_parent_node_->hstart, common_parent_node_->hend);
		}
		return;
	}

	this->ResetSortingIndex(root_);
}

void TreeCommon::RemoveChildNode(CNode* node, bool is_empty_deleted) {
	CBranch* temp_parent = node->parent;
	if (temp_parent == NULL) return;

	int node_index = -1;
	for (int j = 0; j < temp_parent->linked_nodes.size(); ++j)
		if (temp_parent->linked_nodes[j] == node) {
			node_index = j;
			break;
		}
	if (node_index != -1) {
		for (int j = node_index; j < temp_parent->linked_nodes.size() - 1; ++j)
			temp_parent->linked_nodes[j] = temp_parent->linked_nodes[j + 1];
		temp_parent->linked_nodes.resize(temp_parent->linked_nodes.size() - 1);

		if (temp_parent->linked_nodes.size() == 0 && is_empty_deleted)
			this->RemoveChildNode(temp_parent, true);
		else {
			ProgressNodeAndParentData(temp_parent);
			ResetSortingIndex(temp_parent);
		}
	}
}

void TreeCommon::UpdateChildLevel(CBranch* node) {
	std::queue< CNode* > node_queue;
	for (int i = 0; i < node->linked_nodes.size(); ++i)
		node_queue.push(node->linked_nodes[i]);

	while (node_queue.size() > 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();
		if (temp_node->type() == CNode::BRANCH) {
			CBranch* temp_branch = dynamic_cast<CBranch*>(temp_node);
			for (int i = 0; i < temp_branch->linked_nodes.size(); ++i)
				node_queue.push(temp_branch->linked_nodes[i]);
		}
		//temp_node->radius = temp_node->parent->radius / factor_;
		temp_node->set_level(temp_node->parent->level() + 1);
	}

	AssignColor(node, node->hstart, node->hend);
}

int TreeCommon::FindCommonParent(CNode* node, std::vector< int >& node_ids) {
	int value = 0;
	for (int i = 0; i < node_ids.size(); ++i)
		if (node->id() == node_ids[i]) {
			value = 1;
		}

	if (node->type() == CNode::BRANCH) {
		CBranch* branch = dynamic_cast<CBranch*>(node);
		for (int i = 0; i < branch->linked_nodes.size(); ++i)
			value += FindCommonParent(branch->linked_nodes[i], node_ids);
	}
	if (value == node_ids.size()) {
		common_parent_node_ = node;
		value = 0;
	}
	return value;
}

void TreeCommon::ProgressNodeData(CNode* node) {
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
        if (point_index.size() > 1)
            variance[i] = sqrt(variance[i] / (point_index.size() - 1));
        else
		    variance[i] = sqrt(variance[i] / point_index.size());
	node->average_values = average;
	node->variable_variances = variance;
	node->general_variance = 0;
	for (int i = 0; i < variance.size(); ++i)
		node->general_variance += variance[i] * dataset_->var_weights[i];
	node->point_count = point_index.size();
	node->center_pos = center_pos;
}

void TreeCommon::ProgressNodeAndParentData(CNode* node) {
	if (node == NULL || node == common_parent_node_) return;
	ProgressNodeData(node);
	ProgressNodeAndParentData(node->parent);
}

void TreeCommon::GetNodeValues(CNode* node, int var_index, std::vector< float >& values)
{
	std::vector< int > point_index;
	this->Traverse(node, point_index);
	values.resize(point_index.size());
	for (int i = 0; i < point_index.size(); ++i)
		values[i] = dataset_->point_values[point_index[i]][var_index];
}

void TreeCommon::GetConnectionStatus(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status, float& edge_length)
{
    if (dataset_->type() == ScatterPointDataset::POINT_DATA) {
        Utility::VtkTriangulation(nodes, connecting_status, edge_length);
    }
    else {
        connecting_status.resize(nodes.size());
	    for (int i = 0; i < nodes.size(); ++i) {
		    connecting_status[i].resize(nodes.size());
		    connecting_status[i].assign(nodes.size(), false);
	    }

        min_edge_length_ = 1e10;

        for (int i = 0; i < nodes.size() - 1; ++i) {
            std::map< int, int >::iterator iter_one = grid_id_seq_map_.find(nodes[i]->id());
            if (iter_one == grid_id_seq_map_.end()) continue;
            for (int j = i + 1; j < nodes.size(); ++j) {
                std::map< int, int >::iterator iter_two = grid_id_seq_map_.find(nodes[j]->id());
                if (iter_two == grid_id_seq_map_.end()) continue;

                connecting_status[i][j] = node_connecting_status_[iter_one->second][iter_two->second];
                connecting_status[j][i] = connecting_status[i][j];

                if (min_edge_length_ > 1e05)
                    edge_length = sqrt(pow(nodes[i]->center_pos[0] - nodes[j]->center_pos[0], 2) + pow(nodes[i]->center_pos[1] - nodes[j]->center_pos[1], 2));
            }
        }
    }
}
