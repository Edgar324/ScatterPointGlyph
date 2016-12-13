/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef EDGE_NODE_H_
#define EDGE_NODE_H_

#include <vector>
#include <map>
using namespace std;

class MultivariateDataset;

class GaussianNode
{
public:
    GaussianNode() {}
    ~GaussianNode() {}

    vector<double> means;
    vector<double> sigmas;
    vector<double> weights;

    bool fit(vector<double>& val, double alpha = 0.2);
};

class GaussianAggregatedNode
{
public:
    GaussianAggregatedNode() {}
    ~GaussianAggregatedNode() {}

    vector<GaussianNode*> gnodes;

    void CollectIds(MultivariateDataset* dataset, vector<int>& record_ids, double alpha);
};


class EdgeNode
{
public:
    EdgeNode(GaussianAggregatedNode* node_one_t, GaussianAggregatedNode* node_two_t = NULL, EdgeNode* parent_t = NULL);
    ~EdgeNode();

    vector<EdgeNode*> children;

    GaussianAggregatedNode* node_one() { return node_one_; }
    GaussianAggregatedNode* node_two() { return node_two_; }
    int id() { return id_; }

private:
    int id_;
    EdgeNode* parent_;
    GaussianAggregatedNode* node_one_ = NULL;
    GaussianAggregatedNode* node_two_ = NULL;

    static int max_id_;
};

#endif