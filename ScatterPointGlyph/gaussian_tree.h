/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef GAUSSIAN_TREE_H_
#define GAUSSIAN_TREE_H_

#include <vector>
#include <map>
using namespace std;

#include "edge_node.h"

class MultivariateDataset;

class GaussianTree
{
public:
    GaussianTree(MultivariateDataset* dataset_t);
    ~GaussianTree();

    void Build();

private:
    MultivariateDataset* dataset_ = NULL;

    EdgeNode* root_ = NULL;
    map<int, EdgeNode*> id_node_map_;

    void ConstructGaussianNodes(vector<int>& record_ids, vector<GaussianNode*>& gaussian_nodes, vector<vector<bool>>& connecting);
    void Split(vector<GaussianNode*>& gaussian_nodes, vector<vector<bool>>& connecting, vector<vector<int>>& cluster_index);
    void SplitEdgeNode(EdgeNode* node);

    void CollectRecordIds(vector<int>& record_ids, EdgeNode* node, double alpha = 0.2);
};

#endif