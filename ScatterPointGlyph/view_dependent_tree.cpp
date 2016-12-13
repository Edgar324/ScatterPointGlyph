/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "view_dependent_tree.h"
#include "scatter_point_dataset.h"


ViewDependentTree::ViewDependentTree(MultivariateDataset* data)
    : MultiLabelTree(data) {

}

ViewDependentTree::~ViewDependentTree() {

}

void ViewDependentTree::ConstructTree(float left, float right, float bottom, float top, float glyph_radius) {
    //// Step one: reconstruct the root node and separate leaf nodes into two branch nodes.
    //// One contains leaf nodes that are in the viewport, the other contains other leaf nodes.
    //for (int i = 0; i < root_->children.size(); ++i) {
    //    if (root_->children[i]->type() == CNode::BRANCH) {
    //        CBranch* branch = (CBranch*)(root_->children[i]);
    //        branch->children.clear();
    //        id_node_map_.erase(branch->id());
    //        delete branch;
    //    }
    //}
    //root_->children = leaf_nodes;

    //// Step two: split the node that is in the viewport
    //CBranch* view_node = new CBranch;
    //id_node_map_.insert(map<int, CNode*>::value_type(view_node->id(), view_node));
    //CBranch* out_node = new CBranch;
    //id_node_map_.insert(map<int, CNode*>::value_type(out_node->id(), out_node));
    //out_node->left = FLT_MAX;
    //out_node->right = FLT_MAX;
    //out_node->bottom = FLT_MAX;
    //out_node->top = FLT_MAX;
    //out_node->center_pos.resize(2);
    //out_node->center_pos[0] = FLT_MAX;
    //out_node->center_pos[1] = FLT_MAX;

    //for (int i = 0; i < root_->children.size(); ++i) {
    //    float x = root_->children[i]->mean_pos[0] * mv_dataset_->max_pos_range + mv_dataset_->original_pos_ranges[0][0];
    //    float y = root_->children[i]->mean_pos[1] * mv_dataset_->max_pos_range + mv_dataset_->original_pos_ranges[1][0];

    //    if ((x - left) * (x - right) < 0 && (y - bottom) * (y - top) < 0)
    //        view_node->children.push_back(root_->children[i]);
    //    else
    //        out_node->children.push_back(root_->children[i]);
    //}
    //this->ProgressNodeData(view_node);
    //view_node->radius = sqrt((view_node->right - view_node->left) * (view_node->top - view_node->bottom) / EXPECTED_CLUSTER_NUM);

    //root_->children.clear();
    //root_->children.push_back(view_node);
    //root_->children.push_back(out_node);
    //root_->average_dis = 1.0;

    //SplitNodeOnce(view_node->id());
}
