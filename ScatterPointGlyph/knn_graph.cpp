/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "knn_graph.h"
#include <algorithm>
using namespace std;

bool cmp_knn(const KnnNode& lhs, const KnnNode& rhs) {
    return lhs.dis < rhs.dis;
}

KnnGraph::KnnGraph() {

}

KnnGraph::~KnnGraph() {

}

void KnnGraph::BuildGraph(int nearest_neighbors, const vector<vector<double>>& values, const vector<double>& weights, vector<vector<bool>>& connecting) {
    int node_num = values.size();

    vector<vector<double>> dis_mat;
    dis_mat.resize(node_num, vector<double>(node_num, 0));
    for (int i = 0; i + 1 < node_num; i++)
        for (int j = i + 1; j < node_num; j++) {
            double temp_dis = 0;
            for (int k = 0; k < weights.size(); ++k)
                temp_dis += pow(values[i][k] - values[j][k], 2) * weights[k];
            temp_dis = sqrt(temp_dis);
            dis_mat[i][j] = temp_dis;
            dis_mat[j][i] = temp_dis;
        }

    connecting.clear();
    connecting.resize(node_num, vector<bool>(node_num, false));

    for (int i = 0; i < node_num; i++) {
        vector<KnnNode> knn_vec;
        knn_vec.resize(nearest_neighbors + 1);
        for (int j = 0; j < nearest_neighbors + 1; j++)
            knn_vec[j] = KnnNode{ dis_mat[i][j], j };
        make_heap(knn_vec.begin(), knn_vec.end(), cmp_knn);

        for (int j = nearest_neighbors + 1; j < node_num; j++) {
            knn_vec.push_back(KnnNode{ dis_mat[i][j], j });
            push_heap(knn_vec.begin(), knn_vec.end(), cmp_knn);
            pop_heap(knn_vec.begin(), knn_vec.end(), cmp_knn);
            knn_vec.pop_back();
        }

        for (int j = 0; j < nearest_neighbors; j++)
            connecting[i][knn_vec[j + 1].index] = true;
    }
}
