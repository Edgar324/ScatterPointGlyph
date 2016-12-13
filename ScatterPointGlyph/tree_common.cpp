#include "tree_common.h"
#include <iostream>
#include <queue>
#include <time.h>

#include "multivariate_dataset.h"

TreeCommon::TreeCommon(MultivariateDataset* dataset)
    : mv_dataset_(dataset) {
    Initialize();
}

TreeCommon::~TreeCommon() {
    if (root_ != NULL) delete root_;
    id_node_map_.clear();
}

void TreeCommon::Initialize() {
    // Construct root node using leaf_nodes_
    vector<CNode*> temp_nodes(leaf_nodes_.size());
    for (int i = 0; i < leaf_nodes_.size(); i++) 
        temp_nodes[i] = leaf_nodes_[i];
    root_ = new CBranch(mv_dataset_, temp_nodes);
    root_->set_level(0);
    id_node_map_.insert(map<int, CNode*>::value_type(root_->id(), root_));

    for (int i = 0; i < leaf_nodes_.size(); i++) {
        id_node_map_.insert(map<int, CNode*>::value_type(leaf_nodes_[i]->id(), leaf_nodes_[i]));
    }
    root_->Update();

    max_level_ = 1;
}

void TreeCommon::GetNodes(float left, float right, float bottom, float top, vector<CNode*>& nodes) {

    //if (this->type() == GEO_VIEW_DEPENDENT_MULTI_LABEL_TREE || this->type() == KNN_MULTI_LABEL_TREE) {
    //    this->ConstructTree(left, right, bottom, top, glyph_radius);
    //}

    //// Get nodes according to the average distance among all child nodes!!!
    //pre_level_nodes.clear();
    //current_level_nodes.clear();

    //float expected_radius = glyph_radius * 2;
    //float view_center_x = (left + right) / 2;
    //float view_center_y = (top + bottom) / 2;
    //float view_width = right - left;
    //float view_height = top - bottom;

    //if (this->type() == GEO_HIERARCHICAL_TREE) {
    //    pre_level_nodes.push_back(root_);
    //} else if (this->type() == GEO_NCUTS_TREE) {
    //    int temp_level = 0;
    //    vector<CNode*> level_nodes;
    //    do {
    //        temp_level++;
    //        GetNodes(temp_level, level_nodes);
    //    } while (level_nodes.size() < 8);
    //    GetNodes(temp_level - 1, pre_level_nodes);
    //} else {
    //    queue<CNode*> node_queue;
    //    node_queue.push(root_);
    //    while (!node_queue.empty()) {
    //    CNode* node = node_queue.front();
    //    node->is_visible = false;
    //    node_queue.pop();
    //    if (node->type() != CNode::BRANCH) continue;

    //    float center_x = node->center_pos[0] * mv_dataset_->max_pos_range + mv_dataset_->original_pos_ranges[0][0];
    //    float center_y = node->center_pos[1] * mv_dataset_->max_pos_range + mv_dataset_->original_pos_ranges[1][0];
    //    float node_width = (node->right - node->left) * mv_dataset_->max_pos_range;
    //    float node_height = (node->top - node->bottom) * mv_dataset_->max_pos_range;


    //    if (abs(view_center_x - center_x)  < (node_width / 2 + view_width / 2)
    //        && abs(view_center_y - center_y) < (node_height / 2 + view_height / 2)) {
    //        node->is_expanded = true;
    //        if (node->average_dis < expected_radius * 4) {
    //            pre_level_nodes.push_back(node);
    //        } else {
    //            CBranch* branch = (CBranch*)node;
    //            for (int i = 0; i < branch->children.size(); ++i)
    //                node_queue.push(branch->children[i]);
    //        }
    //    }
    //}
    //}

    //for (int i = 0; i < pre_level_nodes.size(); ++i) {
    //    CBranch* branch = (CBranch*)pre_level_nodes[i];
    //    bool is_children_leaf = false;
    //    for (int j = 0; j < branch->children.size(); ++j)
    //        if (branch->children[j]->type() == CNode::LEAF) {
    //            is_children_leaf = true;
    //            break;
    //        }
    //    if (!is_children_leaf) {
    //        for (int j = 0; j < branch->children.size(); ++j) {
    //            CNode* node = branch->children[j];

    //            float center_x = node->center_pos[0] * mv_dataset_->max_pos_range + mv_dataset_->original_pos_ranges[0][0];
    //            float center_y = node->center_pos[1] * mv_dataset_->max_pos_range + mv_dataset_->original_pos_ranges[1][0];
    //            float node_width = (node->right - node->left) * mv_dataset_->max_pos_range;
    //            float node_height = (node->top - node->bottom) * mv_dataset_->max_pos_range;

    //            if (abs(view_center_x - center_x)  < (node_width / 2 + view_width / 2)
    //                && abs(view_center_y - center_y) < (node_height / 2 + view_height / 2))
    //                current_level_nodes.push_back(branch->children[j]);
    //        }
    //    } else {
    //        current_level_nodes.push_back(branch);
    //    }
    //}
    //for (int i = 0; i < current_level_nodes.size(); ++i) {
    //    current_level_nodes[i]->is_expanded = false;
    //    current_level_nodes[i]->is_visible = true;
    //}
}

void TreeCommon::GetNodes(int level, vector<CNode*>& level_nodes) {
    level_nodes.clear();

    queue<CNode*> node_queue;
    node_queue.push(root_);

    while (!node_queue.empty()) {
        CNode* temp_node = node_queue.front();
        node_queue.pop();

        if (temp_node->level() == level && temp_node->type() != CNode::LEAF) {
            level_nodes.push_back(temp_node);
        } else if (temp_node->type() == CNode::BRANCH) {
            CBranch* branch_node = (CBranch*)temp_node;
            bool is_children_all_leaf = true;
            for (int i = 0; i < branch_node->children().size(); ++i)
                if (branch_node->children()[i]->type() == CNode::BRANCH) {
                    is_children_all_leaf = false;
                    break;
                }
            if (is_children_all_leaf) {
                level_nodes.push_back(temp_node);
            } else {
                for (int i = 0; i < branch_node->children().size(); ++i)
                node_queue.push(branch_node->children()[i]);
            }
        }
    }
}

CNode* TreeCommon::GetNode(int node_id) {
    map<int, CNode*>::iterator iter = id_node_map_.find(node_id);
    if (iter != id_node_map_.end()) return iter->second;
    else return NULL;
}

void TreeCommon::run() {
    cout << "Currently not implemented!" << endl;
}

void TreeCommon::Clear() {
    if (root_ == NULL) return;
    /// TODO: implement the clear action
}

//void TreeCommon::GetNodePoints(int node_id, vector<int>& point_ids) {
//    point_ids.clear();
//
//    CNode* node = this->GetNode(node_id);
//    if (node != NULL) this->GetNodePoints(node, point_ids);
//}
//
//void TreeCommon::GetNodePoints(CNode* node, std::vector<int>& point_ids) {
//    point_ids.clear();
//
//	std::queue<CNode*> node_queue;
//	node_queue.push(node);
//
//	while (node_queue.size() != 0) {
//		CNode* temp_node = node_queue.front();
//		node_queue.pop();
//
//		if (temp_node->type() == CNode::BRANCH) {
//			CBranch* branch = dynamic_cast< CBranch* >(temp_node);
//			if (branch != NULL) {
//				for (int i = 0; i < branch->children.size(); ++i)
//					node_queue.push(branch->children[i]);
//			}
//		}
//		else if (temp_node->type() == CNode::LEAF) {
//			CLeaf* leaf = dynamic_cast< CLeaf* >(temp_node);
//			if (leaf != NULL) {
//				for (int i = 0; i < leaf->point_ids.size(); ++i)
//					point_ids.push_back(leaf->point_ids[i]);
//			}
//		}
//	}
//}
//
//void TreeCommon::ConstructRootNode() {
//#ifdef DEBUG_ON
//    cout << "Constructing root node" << endl;
//#endif
//    
//    root_ = new CBranch;
//    root_->set_level(0);
//    root_->radius = 0.5;
//	root_->children.clear();
//    root_->average_dis = 1.0;
//
//    id_node_map_.insert(map<int, CNode*>::value_type(root_->id(), root_));
//
//	float temp_min = 1e10;
//    for (int i = 0; i < mv_dataset_->normalized_point_pos[0].size(); ++i) {
//        if (!mv_dataset_->is_valid[i]) continue;
//
//        CLeaf* temp_leaf = new CLeaf();
//        temp_leaf->point_count = 1;
//        temp_leaf->mean_pos = vector<float>{mv_dataset_->normalized_point_pos[0][i], mv_dataset_->normalized_point_pos[1][i]};
//        temp_leaf->mean_values.resize(mv_dataset_->var_num);
//        for (int j = 0; j < mv_dataset_->var_num; ++j)
//            temp_leaf->mean_values[j] = mv_dataset_->normalized_point_values[j][i];
//        //temp_leaf->mean_values = point_dataset_->normalized_point_values[i];
//        temp_leaf->std_deviations.resize(mv_dataset_->var_num, 0);
//        temp_leaf->parent = root_;
//        temp_leaf->set_level(1);
//        temp_leaf->left = temp_leaf->right = temp_leaf->mean_pos[0];
//        temp_leaf->bottom = temp_leaf->top = temp_leaf->mean_pos[1];
//        temp_leaf->average_dis = FLT_MAX;
//        temp_leaf->point_ids.push_back(i);
//
//        root_->children.push_back(temp_leaf);
//    }
//	this->ProgressNodeData(root_);
//
//#ifdef DEBUG_ON
//    cout << "Root node construction finished!" << endl;
//#endif
//}
//
//void TreeCommon::ConstructLargeScaleRootNode() {
//    #ifdef DEBUG_ON
//    cout << "Constructing root node" << endl;
//#endif
//    
//    root_ = new CBranch;
//    root_->set_level(0);
//    root_->radius = 0.5;
//	root_->children.clear();
//    root_->average_dis = 1.0;
//
//    id_node_map_.insert(map<int, CNode*>::value_type(root_->id(), root_));
//
//    int reso_size = 1000;
//    vector<vector<CLeaf*>> leaf_nodes;
//    leaf_nodes.resize(reso_size);
//    for (int i = 0; i < reso_size; ++i) leaf_nodes[i].resize(reso_size, NULL);
//
//	float temp_min = 1e10;
//    for (int i = 0; i < mv_dataset_->normalized_point_pos[0].size(); ++i) {
//        if (!mv_dataset_->is_valid[i]) continue;
//
//        int x = (int)(mv_dataset_->normalized_point_pos[0][i] * (reso_size - 1));
//        int y = (int)(mv_dataset_->normalized_point_pos[1][i] * (reso_size - 1));
//        CLeaf* temp_leaf = leaf_nodes[y][x];
//        if (temp_leaf == NULL) {
//            temp_leaf = new CLeaf();
//            leaf_nodes[y][x] = temp_leaf;
//            temp_leaf->parent = root_;
//            temp_leaf->set_level(1);
//            root_->children.push_back(temp_leaf);
//        }
//
//        temp_leaf->point_ids.push_back(i);
//    }
//
//    for (int y = 0; y < reso_size; ++y)
//        for (int x = 0; x < reso_size; ++x) 
//            if (leaf_nodes[y][x] != NULL) this->ProgressNodeData(leaf_nodes[y][x]);
//
//	this->ProgressNodeData(root_);
//
//#ifdef DEBUG_ON
//    cout << "Root node construction finished!" << endl;
//#endif
//}
//
//void TreeCommon::ProgressNodeData(CNode* node) {
//	if (node == NULL) return;
//
//	std::vector<int> point_index;
//	this->GetNodePoints(node, point_index);
//
//    if (point_index.size() == 0) {
//        cout << "Error: Node " << node->id() << " have 0 points!" << endl;
//        return;
//    }
//
//    node->point_count = point_index.size();
//	node->mean_values.resize(mv_dataset_->var_num, 0);
//    node->mean_values.assign(mv_dataset_->var_num, 0);
//	node->std_deviations.resize(mv_dataset_->var_num, 0);
//    node->std_deviations.assign(mv_dataset_->var_num, 0);
//	node->mean_pos.resize(2, 0);
//    node->mean_pos.assign(2, 0);
//    node->left = INT_MAX;
//    node->right = INT_MIN;
//    node->top = INT_MIN;
//    node->bottom = INT_MAX;
//
//    float x, y;
//	for (int i = 0; i < point_index.size(); ++i) {
//        x = mv_dataset_->normalized_point_pos[0][point_index[i]];
//        y = mv_dataset_->normalized_point_pos[1][point_index[i]];
//		node->mean_pos[0] += x;
//		node->mean_pos[1] += y;
//        if (x > node->right) node->right = x;
//        if (x < node->left) node->left = x;
//        if (y > node->top) node->top = y;
//        if (y < node->bottom) node->bottom = y;
//
//		for (int j = 0; j < mv_dataset_->var_num; ++j)
//			node->mean_values[j] += mv_dataset_->normalized_point_values[j][point_index[i]];
//	}
//	for (int i = 0; i < mv_dataset_->var_num; ++i)
//		node->mean_values[i] /= point_index.size();
//	node->mean_pos[0] /= point_index.size();
//	node->mean_pos[1] /= point_index.size();
//
//    node->center_pos.resize(2);
//    node->center_pos[0] = (node->left + node->right) / 2;
//    node->center_pos[1] = (node->top + node->bottom) / 2;
//
//	for (int i = 0; i < point_index.size(); ++i)
//		for (int j = 0; j < mv_dataset_->var_num; ++j)
//			node->std_deviations[j] += pow(mv_dataset_->normalized_point_values[j][point_index[i]] - node->mean_values[j], 2);
//	for (int i = 0; i < mv_dataset_->var_num; ++i)
//        if (point_index.size() > 1)
//            node->std_deviations[i] = sqrt(node->std_deviations[i] / (point_index.size() - 1));
//        else
//		    node->std_deviations[i] = sqrt(node->std_deviations[i] / point_index.size());
//}
//
//void TreeCommon::MergeNodes(vector<int>& node_ids) {
//
//}
//
//void TreeCommon::SplitNodeOnce(int node_id) {
//    CNode* node = this->GetNode(node_id);
//    if (node == NULL || node->type() != CNode::BRANCH) {
//        cout << "Error Split Node: " << node_id << endl;
//        return;
//    }
//    CBranch* branch = (CBranch*)node;
//
//	bool is_children_all_leaf = true;
//	for (int i = 0; i < branch->children.size(); ++i)
//		if (branch->children[i]->type() != CNode::LEAF) {
//			is_children_all_leaf = false;
//			break;
//		}
//    if (is_children_all_leaf && branch->children.size() > 2) {
//        SplitNode(branch);
//
//        vector<vector<float>> center_pos;
//        for (int i = 0; i < branch->children.size(); ++i)
//            center_pos.push_back(branch->children[i]->center_pos);
//        branch->average_dis = Utility::GetAverageDistance(center_pos);
//
//        if (branch->level() + 1 > this->max_level_) this->max_level_ = branch->level() + 1;
//
//        for (int i = 0; i < branch->children.size(); ++i) {
//            id_node_map_.insert(map<int, CNode*>::value_type(branch->children[i]->id(), branch->children[i]));
//        }
//    }
//}
//
//void TreeCommon::SplitNodeRecursively(int node_id, float std_dev_threshold) {
//    queue<CNode*> node_queue;
//	node_queue.push(this->GetNode(node_id));
//	while (node_queue.size() > 0) {
//		CNode* temp_node = node_queue.front();
//		node_queue.pop();
//
//		if (temp_node->type() == CNode::BRANCH) {
//			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
//			if (branch == NULL || branch->level() > MAX_ALLOWED_LEVEL) continue;
//
//			float accu_std_dev = 0.0;
//			for (int i = 0; i < mv_dataset_->var_num; ++i) {
//				accu_std_dev += branch->std_deviations[i] * mv_dataset_->var_weights[i];
//			}
//            if ((branch == root_ || accu_std_dev > std_dev_threshold) && branch->point_count > 10) {
//                SplitNode(branch);
//
//                vector<vector<float>> center_pos;
//                for (int i = 0; i < branch->children.size(); ++i)
//                    center_pos.push_back(branch->children[i]->center_pos);
//                branch->average_dis = Utility::GetAverageDistance(center_pos);
//
//                if (branch->level() + 1 > this->max_level_) this->max_level_ = branch->level() + 1;
//
//                for (int i = 0; i < branch->children.size(); ++i) {
//                    id_node_map_.insert(map<int, CNode*>::value_type(branch->children[i]->id(), branch->children[i]));
//                }
//            }
//
//			for (int i = 0; i < branch->children.size(); ++i)
//				node_queue.push(branch->children[i]);
//		}
//	}
//}
//
//
//void TreeCommon::ResetLevel(CNode* node, int level) {
//	if (node == NULL) return;
//	node->set_level(level);
//	if (level > max_level_) max_level_ = level;
//
//	if (node->type() == CNode::BRANCH) {
//		CBranch* branch = dynamic_cast<CBranch*>(node);
//		if (branch != NULL) {
//			for (int i = 0; i < branch->children.size(); ++i)
//				ResetLevel(branch->children[i], level + 1);
//		}
//	}
//}

// ABOVE IS What we really need!

//void TreeCommon::Traverse(int level, std::vector<CNode*>& nodes) {
//    nodes.clear();
//    TraverseLevelNodes(level, root_, nodes);
//	
//	/*if (tree_mode_ == VIEWING_MODE)
//		TraverseLevelNodes(level, root_, nodes);
//	else
//		TraversAllNodes(root_, nodes);*/
//}
//
//void TreeCommon::TraversAllNodes(CNode* root_node, std::vector<CNode*>& nodes) {
//	if (root_node != NULL && root_node->type() == CNode::BRANCH) {
//		CBranch* branch = dynamic_cast<CBranch*>(root_node);
//		branch->is_expanded = true;
//		branch->is_highlighted = false;
//
//		bool is_expandable = true;
//		for (int i = 0; i < branch->linked_nodes.size(); ++i) 
//			if (branch->linked_nodes[i]->type() == CNode::LEAF) {
//				is_expandable = false;
//				break;
//			}
//		if (!is_expandable) {
//			branch->is_expanded = false;
//			nodes.push_back(branch);
//		} else {
//			for (int i = 0; i < branch->linked_nodes.size(); ++i)
//				TraversAllNodes(branch->linked_nodes[branch->sorting_index[i]], nodes);
//		}
//	}
//}
//
//void TreeCommon::TraverseLevelNodes(int level, CNode* root_node, std::vector<CNode*>& nodes) {
//	root_node->is_expanded = true;
//
//	if (root_node->level() == level) {
//		nodes.push_back(root_node);
//		root_node->is_expanded = false;
//	}
//	else {
//		if (root_node->type() == CNode::BRANCH && root_node->level() < level) {
//			CBranch* branch = dynamic_cast<CBranch*>(root_node);
//			bool is_all_leaf = true;
//			for (int i = 0; i < branch->linked_nodes.size(); ++i)
//				if (branch->linked_nodes[i]->type() != CNode::LEAF) {
//					is_all_leaf = false;
//					break;
//				}
//			if (is_all_leaf && level != max_level_) {
//				nodes.push_back(root_node);
//				root_node->is_expanded = false;
//			}
//			else {
//				branch->is_highlighted = false;
//				for (int i = 0; i < branch->linked_nodes.size(); ++i)
//					TraverseLevelNodes(level, branch->linked_nodes[branch->sorting_index[i]], nodes);
//			}
//		} else if (root_node->type() == CNode::LEAF && root_node->level() < level) {
//			nodes.push_back(root_node);
//			root_node->is_expanded = false;
//		}
//	}
//}
//
//void TreeCommon::AssignLeafLevel(CNode* node, int level) {
//	std::queue<CNode*> node_queue;
//	node_queue.push(node);
//
//	while (node_queue.size() != 0) {
//		CNode* temp_node = node_queue.front();
//		node_queue.pop();
//
//		if (temp_node->type() == CNode::BRANCH) {
//			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
//			if (branch != NULL) {
//				for (int i = 0; i < branch->linked_nodes.size(); ++i)
//					node_queue.push(branch->linked_nodes[i]);
//			}
//		}
//		else if (temp_node->type() == CNode::LEAF) {
//			temp_node->set_level(level);
//		}
//	}
//}
//
//void TreeCommon::AssignColor(CNode* node, float hstart, float hend, float factor, bool perm, bool rev) {
//	node->color.setHsl((hstart + hend) / 2 * 255, 0.6 * 255, 0.7 * 255);
//	node->hstart = hstart;
//	node->hend = hend;
//
//	int count = 1;
//	if (node->type() != CNode::LEAF) {
//		CBranch* branch = (CBranch*)node;
//		count = branch->linked_nodes.size();
//	}
//
//	float step = (hend - hstart) / count;
//	std::vector<float> temp_start, temp_end;
//	temp_start.resize(count);
//	temp_end.resize(count);
//
//	for (int i = 0; i < count; ++i){
//		temp_start[i] = hstart + i * step;
//		temp_end[i] = hstart + (i + 1) * step;
//	}
//
//	if (perm) {
//		for (int i = 0; i < count / 2; ++i) {
//			float temp = temp_start[i];
//			temp_start[i] = temp_start[count - i - 1];
//			temp_start[count - i - 1] = temp;
//
//			temp = temp_end[i];
//			temp_end[i] = temp_end[count - i - 1];
//			temp_end[count - i - 1] = temp;
//		}
//	}
//
//	if (rev) {
//		int temp_count = (count - 1) / 4;
//		for (int i = 0; i < temp_count; ++i) {
//			float temp = temp_start[i * 2];
//			temp_start[i * 2] = temp_start[temp_count * 4 - i * 2];
//			temp_start[temp_count * 4 - i * 2] = temp;
//
//			temp = temp_end[i * 2];
//			temp_end[i * 2] = temp_end[temp_count * 4 - i * 2];
//			temp_end[temp_count * 4 - i * 2] = temp;
//		}
//	}
//
//	if (node->type() != CNode::LEAF) {
//		CBranch* branch = (CBranch*)node;
//		for (int i = 0; i < branch->linked_nodes.size(); ++i) {
//			AssignColor(branch->linked_nodes[i], temp_start[i] + (1.0 - factor) / 2.0 * step, temp_end[i] - (1.0 - factor) / 2.0 * step, factor, perm, rev);
//		}
//	}
//}
//
//void TreeCommon::ResetSortingIndex(CNode* node) {
//	std::queue<CNode*> node_queue;
//	node_queue.push(node);
//
//	while (node_queue.size() != 0) {
//		CNode* temp_node = node_queue.front();
//		node_queue.pop();
//
//		if (temp_node->type() == CNode::BRANCH) {
//			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
//			if (branch != NULL) {
//				branch->sorting_index.resize(branch->linked_nodes.size());
//				for (int i = 0; i < branch->sorting_index.size(); ++i)
//					branch->sorting_index[i] = i;
//
//				for (int i = 0; i < branch->linked_nodes.size(); ++i)
//					if (branch->linked_nodes[i]->type() == CNode::BRANCH) node_queue.push(branch->linked_nodes[i]);
//			}
//		}
//	}
//}

//
//void TreeCommon::SortTree(std::vector<int>& node_ids) {
//	int node_count = 0;
//	this->SortNode(root_, node_ids, node_count);
//}
//
//int TreeCommon::SortNode(CNode* node, std::vector<int>& node_ids, int& node_count) {
//	if (node_count > node_ids.size()) return 9999;
//
//	// find all the nodes in the linked_nodes that exist in the node_index
//	if (node->type() == CNode::BRANCH) {
//		CBranch* branch = dynamic_cast<CBranch*>(node);
//		for (int i = 0; i < node_ids.size(); ++i)
//			if (node_ids[i] == node->id()) return i;
//
//		std::vector<bool> is_node_exist;
//		std::vector<int> new_sequence;
//		std::vector<int> node_min_index;
//		node_min_index.resize(branch->linked_nodes.size(), 9999);
//		is_node_exist.resize(branch->linked_nodes.size(), false);
//
//		for (int i = 0; i < branch->linked_nodes.size(); ++i)
//			if (node_count < node_ids.size()) node_min_index[i] = SortNode(branch->linked_nodes[i], node_ids, node_count);
//		for (int i = 0; i < node_min_index.size(); ++i)
//			if (node_min_index[i] < 9999) is_node_exist[i] = true;
//
//		new_sequence.resize(branch->linked_nodes.size());
//		for (int i = 0; i < branch->linked_nodes.size(); ++i) new_sequence[i] = i;
//		Utility::Sort(node_min_index, new_sequence);
//
//		branch->sorting_index.clear();
//		for (int i = 0; i < new_sequence.size(); ++i)
//			if (node_min_index[i] < 9999) 
//				branch->sorting_index.push_back(new_sequence[i]);
//			else 
//				break;
//
//		for (int i = 0; i < is_node_exist.size(); ++i)
//			if (!is_node_exist[i]) branch->sorting_index.push_back(i);
//
//		return node_min_index[0];
//	}
//
//	return 9999;
//}
//
//void TreeCommon::SplitCluster(int cluster_index) {
//	std::map< int, CNode*>::iterator node_iter = id_node_map_.find(cluster_index);
//	if (node_iter == id_node_map_.end()) {
//		std::cout << "Cluster not found!" << std::endl;
//		return;
//	}
//
//	if (node_iter->second->type() != CNode::BRANCH) {
//		std::cout << "Leaf or unknown cluster is not allowed to be split!" << std::endl;
//		return;
//	}
//
//	CBranch* branch = (CBranch*)(node_iter->second);
//	bool is_child_leaf = true;
//	for (int i = 0; i < branch->linked_nodes.size(); ++i)
//		if (branch->linked_nodes[i]->type() != CNode::LEAF) {
//			is_child_leaf = false;
//			break;
//		}
//	if (is_child_leaf) this->SplitNode(branch);
//
//	ResetLevel(branch, branch->level());
//	ResetSortingIndex(branch);
//	AssignColor(branch, branch->hstart, branch->hend);
//
//	/*CBranch* temp_parent = branch->parent;
//	if (temp_parent == NULL) return;
//
//	int node_index = -1;
//	for (int j = 0; j < temp_parent->linked_nodes.size(); ++j)
//		if (temp_parent->linked_nodes[j] == branch) {
//			node_index = j;
//			break;
//		}
//	if (node_index != -1) {
//		int original_size = temp_parent->linked_nodes.size();
//		temp_parent->linked_nodes.resize(temp_parent->linked_nodes.size() + branch->linked_nodes.size() - 1);
//		for (int i = original_size; i > node_index; --i)
//			temp_parent->linked_nodes[i + branch->linked_nodes.size() - 2] = temp_parent->linked_nodes[i];
//		for (int i = 0; i < branch->linked_nodes.size(); ++i) {
//			temp_parent->linked_nodes[i + node_index] = branch->linked_nodes[i];
//			branch->linked_nodes[i]->parent = temp_parent;
//		}
//
//		temp_parent->sorting_index.resize(temp_parent->linked_nodes.size());
//		for (int i = 0; i < temp_parent->sorting_index.size(); ++i)
//			temp_parent->sorting_index[i] = i;
//
//		ResetLevel(temp_parent, temp_parent->level());
//		AssignColor(temp_parent, temp_parent->hstart, temp_parent->hend);
//	}*/
//}
//
//void TreeCommon::MergeClusters(std::vector<int>& cluster_index) {
//	std::vector<CNode*> cluster_nodes;
//	for (int i = 0; i < cluster_index.size(); ++i) {
//		std::map< int, CNode*>::iterator node_iter = id_node_map_.find(cluster_index[i]);
//		if (node_iter == id_node_map_.end()) {
//			std::cout << "Cluster not found!" << std::endl;
//			return;
//		}
//		cluster_nodes.push_back(node_iter->second);
//	}
//
//	bool is_all_branch = true;
//	for (int i = 0; i < cluster_nodes.size(); ++i) {
//		if (cluster_nodes[i]->type() != CNode::BRANCH) {
//			is_all_branch = false;
//			break;
//		}
//	}
//
//	if (is_all_branch) {
//		// find the maximum level node
//		int min_level_index = 10000;
//		CNode* min_level_node = NULL;
//		for (int i = 0; i < cluster_nodes.size(); ++i)
//			if (cluster_nodes[i]->level() < min_level_index) {
//				min_level_index = cluster_nodes[i]->level();
//				min_level_node = cluster_nodes[i];
//			}
//		if (min_level_node != NULL) {
//			CBranch* parent_node = min_level_node->parent;
//			CBranch* new_branch = new CBranch;
//			new_branch->set_level(min_level_node->level());
//			new_branch->radius = min_level_node->radius;
//			new_branch->parent = parent_node;
//			parent_node->linked_nodes.push_back(new_branch);
//			id_node_map_.insert(std::map< int, CNode*>::value_type(new_branch->id(), new_branch));
//
//			this->FindCommonParent(root_, cluster_index);
//
//			for (int i = 0; i < cluster_nodes.size(); ++i) {
//				RemoveChildNode(cluster_nodes[i], true);
//				CBranch* branch = dynamic_cast<CBranch*>(cluster_nodes[i]);
//				for (int j = 0; j < branch->linked_nodes.size(); ++j)
//					new_branch->linked_nodes.push_back(branch->linked_nodes[j]);
//			}
//
//			ResetSortingIndex(new_branch);
//			UpdateChildLevel(new_branch->parent);
//			ProgressNodeAndParentData(new_branch);
//			AssignColor(common_parent_node_, common_parent_node_->hstart, common_parent_node_->hend);
//		}
//		return;
//	}
//
//	this->ResetSortingIndex(root_);
//}
//
//void TreeCommon::RemoveChildNode(CNode* node, bool is_empty_deleted) {
//	CBranch* temp_parent = node->parent;
//	if (temp_parent == NULL) return;
//
//	int node_index = -1;
//	for (int j = 0; j < temp_parent->linked_nodes.size(); ++j)
//		if (temp_parent->linked_nodes[j] == node) {
//			node_index = j;
//			break;
//		}
//	if (node_index != -1) {
//		for (int j = node_index; j < temp_parent->linked_nodes.size() - 1; ++j)
//			temp_parent->linked_nodes[j] = temp_parent->linked_nodes[j + 1];
//		temp_parent->linked_nodes.resize(temp_parent->linked_nodes.size() - 1);
//
//		if (temp_parent->linked_nodes.size() == 0 && is_empty_deleted)
//			this->RemoveChildNode(temp_parent, true);
//		else {
//			ProgressNodeAndParentData(temp_parent);
//			ResetSortingIndex(temp_parent);
//		}
//	}
//}
//
//void TreeCommon::UpdateChildLevel(CBranch* node) {
//	std::queue<CNode*> node_queue;
//	for (int i = 0; i < node->linked_nodes.size(); ++i)
//		node_queue.push(node->linked_nodes[i]);
//
//	while (node_queue.size() > 0) {
//		CNode* temp_node = node_queue.front();
//		node_queue.pop();
//		if (temp_node->type() == CNode::BRANCH) {
//			CBranch* temp_branch = dynamic_cast<CBranch*>(temp_node);
//			for (int i = 0; i < temp_branch->linked_nodes.size(); ++i)
//				node_queue.push(temp_branch->linked_nodes[i]);
//		}
//		//temp_node->radius = temp_node->parent->radius / factor_;
//		temp_node->set_level(temp_node->parent->level() + 1);
//	}
//
//	AssignColor(node, node->hstart, node->hend);
//}
//
//int TreeCommon::FindCommonParent(CNode* node, std::vector<int>& node_ids) {
//	int value = 0;
//	for (int i = 0; i < node_ids.size(); ++i)
//		if (node->id() == node_ids[i]) {
//			value = 1;
//		}
//
//	if (node->type() == CNode::BRANCH) {
//		CBranch* branch = dynamic_cast<CBranch*>(node);
//		for (int i = 0; i < branch->linked_nodes.size(); ++i)
//			value += FindCommonParent(branch->linked_nodes[i], node_ids);
//	}
//	if (value == node_ids.size()) {
//		common_parent_node_ = node;
//		value = 0;
//	}
//	return value;
//}
//
//void TreeCommon::GetNodeValues(CNode* node, int var_index, std::vector<float>& values)
//{
//	std::vector<int> point_index;
//	this->GetNodePoints(node, point_index);
//	values.resize(point_index.size());
//	for (int i = 0; i < point_index.size(); ++i)
//		values[i] = point_dataset_->normalized_point_values[point_index[i]][var_index];
//}
//
//void TreeCommon::GetConnectionStatus(std::vector<CNode*>& nodes, std::vector<std::vector<bool>>& connecting_status, float& edge_length)
//{
//    if (point_dataset_->type() == ScatterPointDataset::POINT_DATA) {
//        Utility::VtkTriangulation(nodes, connecting_status, edge_length);
//    }
//    else {
//        connecting_status.resize(nodes.size());
//	    for (int i = 0; i < nodes.size(); ++i) {
//		    connecting_status[i].resize(nodes.size());
//		    connecting_status[i].assign(nodes.size(), false);
//	    }
//
//        min_edge_length_ = 1e10;
//
//        for (int i = 0; i < nodes.size() - 1; ++i) {
//            std::map< int, int >::iterator iter_one = grid_id_seq_map_.find(nodes[i]->id());
//            if (iter_one == grid_id_seq_map_.end()) continue;
//            for (int j = i + 1; j < nodes.size(); ++j) {
//                std::map< int, int >::iterator iter_two = grid_id_seq_map_.find(nodes[j]->id());
//                if (iter_two == grid_id_seq_map_.end()) continue;
//
//                connecting_status[i][j] = node_connecting_status_[iter_one->second][iter_two->second];
//                connecting_status[j][i] = connecting_status[i][j];
//
//                if (min_edge_length_ > 1e05)
//                    edge_length = sqrt(pow(nodes[i]->mean_pos[0] - nodes[j]->mean_pos[0], 2) + pow(nodes[i]->mean_pos[1] - nodes[j]->mean_pos[1], 2));
//            }
//        }
//    }
//}
