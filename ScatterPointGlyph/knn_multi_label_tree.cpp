/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "knn_multi_label_tree.h"
#include "cnode.h"
#include "multivariate_dataset.h"
#include "knn_graph.h"
#include "multi_label_processor.h"

KnnMultiLabelTree::KnnMultiLabelTree(MultivariateDataset* data, int nearest_neighbors /*= 5*/) 
    : TreeCommon(data) {
    processor_ = new MultiLabelProcessor;
}

KnnMultiLabelTree::~KnnMultiLabelTree() {

}

void KnnMultiLabelTree::SplitNode(CBranch* root_node) {
    const vector<CNode*>& children = root_node->children();
    vector<vector<double>> values;
    values.resize(children.size());
    for (int i = 0; i < children.size(); i++)
        values[i] = children[i]->mean_values();

    vector<vector<bool>> connecting;
    KnnGraph knn_solver;
    knn_solver.BuildGraph(nearest_neighbors_, values, mv_dataset_->var_weights(), connecting);

    vector<vector<int>> clusters;
    double radius = 0.3;
    processor_->GenerateCluster(values, mv_dataset_->var_weights(), connecting, radius, clusters);

    vector<CNode*> new_children(clusters.size());
    for (int i = 0; i < clusters.size(); i++) {
        vector<CNode*> temp_nodes(clusters[i].size());
        for (int j = 0; j < clusters[i].size(); j++)
            temp_nodes[j] = children[clusters[i][j]];
        new_children[i] = new CBranch(mv_dataset_, temp_nodes);
        new_children[i]->Update(true, false);
    }
    root_node->ClearChildren();
    root_node->AppendChildren(new_children);
}

void KnnMultiLabelTree::AutoConstructTree(float std_dev_threshold) {
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
    root_->set_level(0);
    SplitNode(root_);
}

void KnnMultiLabelTree::ConstructTree(float left, float right, float bottom, float top) {

}

void KnnMultiLabelTree::Clear() {

}

void KnnMultiLabelTree::BuildLeafs() {
    BuildOnRecords();
}

void KnnMultiLabelTree::BuildOnKwayPartitions() {

}

void KnnMultiLabelTree::BuildOnRecords() {
    leaf_nodes_.resize(mv_dataset_->record_num());
    for (int i = 0; i < mv_dataset_->record_num(); ++i) {
        leaf_nodes_[i] = new CLeaf(mv_dataset_, i);
        leaf_nodes_[i]->Update(true, false);

        id_node_map_.insert(map<int, CNode*>::value_type(leaf_nodes_[i]->id(), leaf_nodes_[i]));
    }
}

