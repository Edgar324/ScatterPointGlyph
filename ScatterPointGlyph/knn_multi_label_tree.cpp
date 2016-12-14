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

void KnnMultiLabelTree::SplitNode(CBranch* node) {
    const vector<CNode*>& children = node->children();
    vector<vector<double>> values;
    values.resize(children.size());
    for (int i = 0; i < children.size(); i++)
        values[i] = children[i]->mean_values();

    vector<vector<bool>> connecting;
    KnnGraph knn_solver;
    knn_solver.BuildGraph(nearest_neighbors_, values, mv_dataset_->var_weights(), connecting);

    vector<vector<int>> clusters;
    double radius = 0.2;
    processor_->GenerateCluster(values, mv_dataset_->var_weights(), connecting, radius, clusters);

}

void KnnMultiLabelTree::AutoConstructTree(float std_dev_threshold) {

}

void KnnMultiLabelTree::ConstructTree(float left, float right, float bottom, float top) {

}

void KnnMultiLabelTree::Clear() {

}

void KnnMultiLabelTree::BuildLeafs() {

}

void KnnMultiLabelTree::BuildOnKwayPartitions() {

}

void KnnMultiLabelTree::BuildOnRecords() {

}

