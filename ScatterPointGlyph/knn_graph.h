/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef KNN_GRAPH_H_
#define KNN_GRAPH_H_

#include <vector>
using namespace std;

class MultivariateDataset;

struct KnnNode {
    double dis;
    int index;

    KnnNode& operator=(const KnnNode& knn) {
        dis = knn.dis;
        index = knn.index;

        return *this;
    }
};

class KnnGraph
{
public:
    KnnGraph();
    ~KnnGraph();

    void BuildGraph(int nearest_neighbors, const vector<vector<double>>& values, const vector<double>& weights, vector<vector<bool>>& connecting);
    double max_distance() { return max_distance_; }

private:
    double max_distance_ = 0;

};

#endif