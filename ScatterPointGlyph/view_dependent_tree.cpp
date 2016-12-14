/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "view_dependent_tree.h"
#include "scatter_point_dataset.h"
#include "multivariate_dataset.h"


ViewDependentTree::ViewDependentTree(MultivariateDataset* data)
    : MultiLabelTree(data) {

}

ViewDependentTree::~ViewDependentTree() {

}

void ViewDependentTree::ConstructTree(float left, float right, float bottom, float top) {
    left_ = left;
    right_ = right;
    bottom_ = bottom;
    top_ = top;

    // Step 1: Clear all the nodes
    if (root_ != NULL) {
        delete root_;
        root_ = NULL;
        leaf_nodes_.clear();
        id_node_map_.clear();
    }

    // Step 2: Build leaf nodes that are in the viewport
    BuildLeafs();

    // Step 3: split the leaf nodes
    root_ = new CBranch(mv_dataset_, leaf_nodes_);
    SplitNode(root_);
}

void ViewDependentTree::BuildLeafs() {
    int node_count = 0;
    for (int i = 0; i < mv_dataset_->record_num(); i++) {
        double* pos = mv_dataset_->GetRecordPos(i);
        if ((pos[0] - left_) * (pos[0] - right_) <= 0
            && (pos[1] - bottom_) * (pos[1] - top_) <= 0)
            node_count++;
    }

    if (node_count < SCALE_NODE_SIZE) {
        BuildOnRecords();
        is_pixel_used_ = false;
    } else {
        BuildOnPixels();
        is_pixel_used_ = true;
    }
}

void ViewDependentTree::BuildOnPixels() {
    int imw = LEAF_RESO_SIZE;
    int imh = (top_ - bottom_) / (right_ - left_) * imw;

    vector<CLeaf*> grid_nodes;
    grid_nodes.resize(imw * imh, NULL);

    for (int i = 0; i < mv_dataset_->record_num(); ++i) {
        double* record_pos = mv_dataset_->GetRecordPos(i);

        if (!((record_pos[0] - left_) * (record_pos[0] - right_) <= 0
            && (record_pos[1] - bottom_) * (record_pos[1] - top_) <= 0)) continue;

        int x = (int)(record_pos[0] * (imw - 1));
        int y = (int)(record_pos[1] * (imh - 1));
        int grid_index = y * imw + x;
        
        if (grid_nodes[grid_index] == NULL) {
            grid_nodes[grid_index] = new CLeaf(mv_dataset_, i);

            id_node_map_.insert(map<int, CNode*>::value_type(grid_nodes[grid_index]->id(), grid_nodes[grid_index]));
        } else {
            grid_nodes[grid_index]->AppendRecord(i);
        }
    }
    
    leaf_nodes_.clear();
    for (int i = 0; i < grid_nodes.size(); i++)
        if (grid_nodes[i] != NULL) {
            grid_nodes[i]->Update();
            leaf_nodes_.push_back(grid_nodes[i]);
        }
}

void ViewDependentTree::BuildOnRecords() {
    leaf_nodes_.clear();
    for (int i = 0; i < mv_dataset_->record_num(); ++i) {
        double* record_pos = mv_dataset_->GetRecordPos(i);

        if (!((record_pos[0] - left_) * (record_pos[0] - right_) <= 0
            && (record_pos[1] - bottom_) * (record_pos[1] - top_) <= 0)) continue;

        CLeaf* leaf = new CLeaf(mv_dataset_, i);
        leaf->Update();
        leaf_nodes_.push_back(leaf);

        id_node_map_.insert(map<int, CNode*>::value_type(leaf->id(), leaf));
    }
}
