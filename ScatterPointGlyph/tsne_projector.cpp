/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "tsne_projector.h"
#include "tsne.h"

TsneProjector::TsneProjector() {

}

TsneProjector::~TsneProjector() {

}

void TsneProjector::Project(vector<vector<float>>& values, vector<vector<float>>& proj_values) {
    int N = values[0].size();
    int D = values.size();
    int no_dims = 2;
    double perplexity = 30;
    if (perplexity > values[0].size() / 5) perplexity = values[0].size() / 5;
    double theta = 0.5;
    int rand_seed = 12345;
    bool skip_random_init = false;

    proj_values.resize(no_dims);
    for (int i = 0; i < no_dims; ++i)
        proj_values[i].resize(N);

    vector<double> X, Y;
    X.resize(N * D);
    Y.resize(no_dims * N);

    int accu_count = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < D; ++j) {
            X[accu_count] = values[j][i];
            accu_count++;
        }
    }

    TSNE sne;
    sne.run(X.data(), N, D, Y.data(), no_dims, perplexity, theta, rand_seed, skip_random_init);

    accu_count = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < no_dims; ++j) {
            proj_values[j][i] = Y[accu_count] * 1000;
            accu_count++;
        }
    }
}