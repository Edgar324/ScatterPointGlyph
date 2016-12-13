/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "gaussian_tree.h"
#include "multivariate_dataset.h"
#include "gmm.h"

GaussianTree::GaussianTree(MultivariateDataset* dataset_t) 
    : dataset_(dataset_t) {
    
}

GaussianTree::~GaussianTree() {

}

void GaussianTree::Build() {
    //// Step 1: Construct root node
    //vector<int> record_ids;
    //record_ids.resize(dataset_->record_num);
    //for (int i = 0; i < dataset_->record_num; i++) record_ids[i] = i;
    //vector<GaussianNode*> gaussian_nodes;
    //vector<vector<bool>> connecting;
    //ConstructGaussianNodes(record_ids, gaussian_nodes, connecting);

    //GaussianAggregatedNode* agg_node = new GaussianAggregatedNode;
    //agg_node->gnodes = gaussian_nodes;
    //root_ = new EdgeNode(agg_node);
    //id_node_map_.insert(map<int, EdgeNode*>::value_type(root_->id(), root_));

    //// Step 2: Construct the first level representation
    //SplitEdgeNode(root_);
}

void GaussianTree::ConstructGaussianNodes(vector<int>& record_ids, vector<GaussianNode*>& gaussian_nodes, vector<vector<bool>>& connecting) {
    //const int gaussians = 1;
    //const size_t maxIterations = 250;
    //const double tolerance = 1e-10;
    //vector<double> weight;
    //vector<double> mu;
    //vector<double> sigma;
    //weight.resize(gaussians);
    //mu.resize(gaussians);
    //sigma.resize(gaussians);

    //// Step 1: Set up GMM models for each variable
    //vector<vector<double>> gmm_models;
    //gmm_models.resize(dataset_->var_num);
    //for (int i = 0; i < dataset_->var_num; i++)
    //    gmm_models[i].resize(gaussians * 3);
    //vector<double> var_values;
    //var_values.resize(record_ids.size());
    //for (int i = 0; i < dataset_->var_num; i++) {
    //    double min_value = 1e10;
    //    double max_value = -1e10;
    //    for (int j = 0; j < record_ids.size(); j++) {
    //        float temp_value = dataset_->var_values[i][record_ids[j]];
    //        var_values[j] = temp_value;
    //        if (temp_value > max_value) max_value = temp_value;
    //        if (temp_value < min_value) min_value = temp_value;
    //    }

    //    for (int j = 0; j < gaussians; j++) {
    //        weight[j] = 1.0 / gaussians;
    //        mu[j] = (max_value - min_value) / gaussians * (j + 0.5) + min_value;
    //        sigma[j] = (max_value - min_value) / (gaussians * 2);
    //    }
    //    GMM gmm(gaussians, weight.data(), mu.data(), sigma.data(), maxIterations, tolerance);
    //    gmm.estimate(var_values.data(), var_values.size());

    //    for (int j = 0; j < gaussians; j++) {
    //        gmm_models[i][j * 3] = gmm.getMean(j);
    //        gmm_models[i][j * 3 + 1] = gmm.getVar(j);
    //        gmm_models[i][j * 3 + 2] = gmm.getMixCoefficient(j);
    //    }
    //}

    //// Step 2: Create Gaussian node based on GMM models
    //// Step 3: Create connection map
    //int gaussian_node_size = 1;
    //for (int i = 0; i < gmm_models.size(); i++)
    //    gaussian_node_size *= (gmm_models[i].size() / 3);
    //connecting.resize(gaussian_node_size, vector<bool>(gaussian_node_size, false));

    //vector<int> accu_model_size = vector<int>(gmm_models.size(), 0);
    //vector<int> mul_model_size = vector<int>(gmm_models.size(), 0);
    //int temp_mul = 1;
    //for (int i = 0; i < gmm_models.size(); i++) {
    //    mul_model_size[i] = temp_mul;
    //    temp_mul *= gmm_models[i].size() / 3;
    //}

    //int var_num = gmm_models.size();
    //int current_size = 0;
    //while (true) {
    //    accu_model_size[0]++;
    //    int temp_index = 0;
    //    while (temp_index + 1 < var_num && accu_model_size[temp_index] >= gmm_models[temp_index].size() / 3) {
    //        accu_model_size[temp_index] = 0;
    //        accu_model_size[temp_index + 1]++;
    //        temp_index++;
    //    }
    //    if (accu_model_size[var_num - 1] < gmm_models[var_num - 1].size() / 3) break;

    //    GaussianNode* gnode = new GaussianNode;
    //    for (int i = 0; i < var_num; i++) {
    //        gnode->means.push_back(gmm_models[i][accu_model_size[i] * 3]);
    //        gnode->sigmas.push_back(gmm_models[i][accu_model_size[i] * 3 + 1]);
    //        gnode->weights.push_back(gmm_models[i][accu_model_size[i] * 3 + 2]);
    //    }
    //    gaussian_nodes.push_back(gnode);

    //    for (int i = 0; i < var_num; i++) {
    //        if (accu_model_size[i] + 1 < gmm_models[i].size() / 3) {
    //            int next_index = current_size + mul_model_size[i];
    //            connecting[current_size][next_index] = true;
    //            connecting[next_index][current_size] = true;
    //        }
    //    }

    //    current_size++;
    //}
}

void GaussianTree::Split(vector<GaussianNode*>& gaussian_nodes, vector<vector<bool>>& connecting, vector<vector<int>>& cluster_index) {
    // Split nodes using multi-label optimization
    /*GraphCluster cluster;
    cluster.SetData(gaussian_nodes, connecting);
    double max_range = -1e10, min_range = 1e10;
    for (int i = 0; i + 1 < gaussian_nodes.size(); i++)
        for (int j = i + 1; j < gaussian_nodes.size(); j++)
            for (int k = 0; k < gaussian_nodes[i]->means.size(); k++) {
                double temp_range = abs(gaussian_nodes[i]->means[k] - gaussian_nodes[j]->means[k]);
                if (temp_range > max_range) max_range = temp_range;
                if (temp_range < min_range) min_range = temp_range;
            }
    float child_radius = (min_range + max_range) / 4;
    cluster.SetLabelEstimationRadius(child_radius);
    cluster.GenerateCluster();

    vector<int>& node_labels = cluster.GetResultLabel();*/
}

void GaussianTree::SplitEdgeNode(EdgeNode* node) {
    //// Step 1: Collect Gaussian nodes under the edge
    //vector<int> record_ids;
    //if (node != NULL) CollectRecordIds(record_ids, node);
    //vector<GaussianNode*> gaussian_nodes;
    //vector<vector<bool>> connecting;
    //ConstructGaussianNodes(record_ids, gaussian_nodes, connecting);

    //// Step 2: Split the nodes based on multi-label optimization
    //vector<vector<int>> cluster_index;
    //Split(gaussian_nodes, connecting, cluster_index);
    //vector<GaussianAggregatedNode*> agg_nodes;
    //for (int i = 0; i < cluster_index.size(); i++) {
    //    GaussianAggregatedNode* temp_node = new GaussianAggregatedNode;
    //    for (int j = 0; j < cluster_index[i].size(); j++)
    //        temp_node->gnodes.push_back(gaussian_nodes[cluster_index[i][j]]);
    //    agg_nodes.push_back(temp_node);
    //}

    //// Step 3: Create new aggregated Gaussian nodes and new edges
    //for (int i = 0; i < cluster_index.size(); i++)
    //    for (int j = i + 1; j < cluster_index.size(); j++) {
    //        GaussianAggregatedNode* node_one = agg_nodes[i];
    //        GaussianAggregatedNode* node_two = agg_nodes[j];
    //        
    //        // Check connecting
    //        bool is_connecting = false;
    //        for (int ci = 0; ci < cluster_index[i].size(); ci++)
    //            for (int cj = 0; cj < cluster_index[j].size(); cj++)
    //                if (connecting[cluster_index[i][ci]][cluster_index[j][cj]]) {
    //                    is_connecting = true;
    //                    break;
    //                }

    //        if (!is_connecting) continue;

    //        EdgeNode* temp_edge = new EdgeNode(node_one, node_two, node);
    //        node->children.push_back(temp_edge);
    //    }
}

void GaussianTree::CollectRecordIds(vector<int>& record_ids, EdgeNode* node, double alpha) {
    /*record_ids.clear();

    vector<int> ids_one, ids_two;
    if (node->node_one() != NULL) node->node_one()->CollectIds(dataset_, ids_one, alpha);
    if (node->node_two() != NULL) node->node_two()->CollectIds(dataset_, ids_two, alpha);

    sort(ids_one.begin(), ids_one.end());
    sort(ids_two.begin(), ids_two.end());

    int index_one = 0, index_two = 0;
    int current_id = -1;
    while (true) {
        while (index_one < ids_one.size() && ids_one[index_one] <= current_id) index_one++;
        while (index_two < ids_two.size() && ids_two[index_two] <= current_id) index_two++;
        if (index_one >= ids_one.size() || index_two >= ids_two.size()) break;
        if (ids_one[index_one] <= ids_two[index_two]) {
            record_ids.push_back(ids_one[index_one]);
            current_id = ids_one[index_one];
            index_one++;
        } else {
            record_ids.push_back(ids_two[index_two]);
            current_id = ids_two[index_two];
            index_two++;
        }
    }

    while (index_one < ids_one.size()) {
        record_ids.push_back(ids_one[index_one]);
        index_one++;
    }

    while (index_two < ids_two.size()) {
        record_ids.push_back(ids_two[index_two]);
        index_two++;
    }*/
}
