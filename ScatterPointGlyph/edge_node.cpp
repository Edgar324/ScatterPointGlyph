/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "edge_node.h"
#include "multivariate_dataset.h"

#define SQRT_2PI 2.5066282746310005024147107274575

int EdgeNode::max_id_ = 0;

bool GaussianNode::fit(vector<double>& val, double alpha) {
    bool is_fit = true;
    for (int i = 0; i < val.size(); i++) {
        double gprob = 1.0 / (sigmas[i] * SQRT_2PI) * exp(-1 * pow(val[i] - means[i], 2) / (2 * sigmas[i] * sigmas[i]));
        if (gprob < alpha) {
            return false;
        }
    }
    return true;
}

EdgeNode::EdgeNode(GaussianAggregatedNode* node_one_t, GaussianAggregatedNode* node_two_t, EdgeNode* parent_t) 
    : node_one_(node_one_t), node_two_(node_two_t), parent_(parent_t) {
    this->id_ = max_id_++;
}

EdgeNode::~EdgeNode() {

}

void GaussianAggregatedNode::CollectIds(MultivariateDataset* dataset, vector<int>& record_ids, double alpha) {
    /*vector<double> val = vector<double>(dataset->var_num());

    bool is_fit = false;
    for (int i = 0; i < dataset->record_num(); i++) {
        is_fit = false;
        for (int j = 0; j < dataset->var_num(); j++)
            val[j] = dataset->var_values()[j][i];
        for (int j = 0; j < gnodes.size(); j++)
            if (gnodes[j]->fit(val, alpha)) {
                is_fit = true;
                break;
            }
        if (is_fit) record_ids.push_back(i);
    }*/
}
