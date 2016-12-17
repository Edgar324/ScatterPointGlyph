/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "knn_multi_label_tree.h"
#include <iostream>
#include <algorithm>
using namespace std;
#include <metis.h>
#include "cnode.h"
#include "multivariate_dataset.h"
#include "knn_graph.h"
#include "multi_label_processor.h"

bool cmp_knn(const KnnNode& lhs, const KnnNode& rhs) {
    return lhs.dis < rhs.dis;
}

KnnMultiLabelTree::KnnMultiLabelTree(MultivariateDataset* data, int nearest_neighbors /*= 5*/) 
    : TreeCommon(data) {
    processor_ = new MultiLabelProcessor;

    // Initialize the root tree
    // Step 1: Build leaf nodes that are in the viewport
    current_records_.resize(mv_dataset_->record_num());
    for (int i = 0; i < current_records_.size(); i++)
        current_records_[i] = i;
    BuildLeafs();
    // Step 2: split the leaf nodes
    root_ = new CBranch(mv_dataset_, leaf_nodes_);
    exploration_path_.push_back(vector<CNode*>(1, root_));
    SplitNode(root_);
    exploration_path_.push_back(root_->children());
}

KnnMultiLabelTree::~KnnMultiLabelTree() {

}

void KnnMultiLabelTree::GetNodes(int level, vector<CNode*>& level_nodes) {
    level_nodes.clear();

    if (level >= exploration_path_.size()) {
        cout << "Trying to get nodes which have a higher level than the exploration path!" << endl;
        return;
    }

    level_nodes = exploration_path_[level];
}

void KnnMultiLabelTree::GetNodes(float left, float right, float bottom, float top, vector<CNode*>& nodes) {
    cout << "View dependent GetNodes is not supported in knn multilabel tree" << endl;
    return;
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
    double model_radius = knn_solver.max_distance() * radius_scale_;
    processor_->GenerateCluster(values, mv_dataset_->var_weights(), connecting, model_radius, clusters);

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

void KnnMultiLabelTree::SplitNodes(int level, vector<CBranch*>& nodes) {
    // Step 1: Clear all nodes after level
    for (int i = level + 1; i < exploration_path_.size(); i++) {
        for (int j = 0; j < exploration_path_[i].size(); j++)
            delete exploration_path_[i][j];
        exploration_path_[i].clear();
    }
    exploration_path_.resize(level + 1);

    // Step 2: Collect all leaf nodes and records
    vector<CLeaf*> temp_leaf_ndoes;
    for (int i = 0; i < nodes.size(); i++) {
        const vector<CNode*> children = nodes[i]->children();
        for (int j = 0; j < children.size(); j++)
            temp_leaf_ndoes.push_back((CLeaf*)children[j]);
    }
    current_records_.clear();
    for (int i = 0; i < temp_leaf_ndoes.size(); i++) {
        const vector<int>& record_ids = temp_leaf_ndoes[i]->record_ids();
        int current_size = current_records_.size();
        current_records_.resize(current_size + record_ids.size());
        for (int j = 0; j < record_ids.size(); j++)
            current_records_[current_size + j] = record_ids[j];
    }

    // Step 3: Generate level nodes
    BuildLeafs();

    vector<vector<double>> values;
    values.resize(leaf_nodes_.size());
    for (int i = 0; i < leaf_nodes_.size(); i++)
        values[i] = leaf_nodes_[i]->mean_values();

    vector<vector<bool>> connecting;
    KnnGraph knn_solver;
    knn_solver.BuildGraph(nearest_neighbors_, values, mv_dataset_->var_weights(), connecting);

    vector<vector<int>> clusters;
    double model_radius = knn_solver.max_distance() * radius_scale_;
    processor_->GenerateCluster(values, mv_dataset_->var_weights(), connecting, model_radius, clusters);

    vector<CNode*> new_children(clusters.size());
    for (int i = 0; i < clusters.size(); i++) {
        vector<CNode*> temp_nodes(clusters[i].size());
        for (int j = 0; j < clusters[i].size(); j++)
            temp_nodes[j] = leaf_nodes_[clusters[i][j]];
        new_children[i] = new CBranch(mv_dataset_, temp_nodes);
        new_children[i]->Update(true, false);
    }
    exploration_path_.push_back(new_children);
}

void KnnMultiLabelTree::AutoConstructTree(float std_dev_threshold) {
    cout << "KnnMultiLabelTree : AutoConstructTree currently not supported!" << endl;
    return;
}

void KnnMultiLabelTree::ConstructTree(float left, float right, float bottom, float top) {
    
}

void KnnMultiLabelTree::Clear() {

}

void KnnMultiLabelTree::BuildLeafs() {
    if (current_records_.size() < 2000)
        BuildOnRecords();
    else
        BuildOnKwayPartitions();
}

void KnnMultiLabelTree::BuildOnKwayPartitions() {
    // Step 1: Build KNN graph
    // Step 2: K-way partitioning
    idx_t nvtxs = current_records_.size();
    idx_t ncon = 1;
    idx_t nparts = 1000;
    idx_t objval;
    vector<idx_t> xajd, adjncy, adjwgt, part;

    vector<idx_t> edges;
    vector<double> edge_distance;

    vector<double> distance(nvtxs);
    vector<KnnNode> knn_vec(nearest_neighbors_ + 1);
    for (int i = 0; i < nvtxs; i++) {
        for (int j = 0; j < nvtxs; j++)
            distance[j] = mv_dataset_->Distance(current_records_[i], current_records_[j]);

        for (int j = 0; j < nearest_neighbors_ + 1; j++)
            knn_vec[j] = KnnNode{ distance[j], j };
        make_heap(knn_vec.begin(), knn_vec.end(), cmp_knn);

        for (int j = nearest_neighbors_ + 1; j < nvtxs; j++) {
            knn_vec.push_back(KnnNode{ distance[j], j });
            push_heap(knn_vec.begin(), knn_vec.end(), cmp_knn);
            pop_heap(knn_vec.begin(), knn_vec.end(), cmp_knn);
            knn_vec.pop_back();
        }

        for (int j = 0; j < nearest_neighbors_ + 1; j++) {
            if (knn_vec[j].index == i) continue;
            edges.push_back(i);
            edges.push_back(knn_vec[j].index);
            edge_distance.push_back(knn_vec[j].dis);
        }
    }

    for (int i = 0; i < nvtxs; i++) {
        xajd.push_back(adjncy.size());
        for (int j = 0; j < edges.size(); j += 2) {
            if (edges[j] == i) {
                adjncy.push_back(edges[j + 1]);
                adjwgt.push_back((idx_t)(edge_distance[j / 2] * 10000));
            } else if (edges[j + 1] == i) {
                adjncy.push_back(edges[j]);
                adjwgt.push_back((idx_t)(edge_distance[j / 2] * 10000));
            }
        }
    }
    xajd.push_back(adjncy.size());

    part.resize(nvtxs);
    METIS_PartGraphRecursive(&nvtxs, &ncon, xajd.data(), adjncy.data(), NULL, NULL, adjwgt.data(), &nparts, NULL, NULL, NULL, &objval, part.data());
    cout << "KWay finished" << endl;

    // Step 3: Build leaf nodes
    vector<vector<int>> leaf_records(nparts);
    for (int i = 0; i < part.size(); i++)
        leaf_records[part[i]].push_back(current_records_[i]);

    leaf_nodes_.resize(nparts);
    for (int i = 0; i < nparts; i++) {
        leaf_nodes_[i] = new CLeaf(mv_dataset_, leaf_records[i]);
        leaf_nodes_[i]->Update(true, false);

        id_node_map_.insert(map<int, CNode*>::value_type(leaf_nodes_[i]->id(), leaf_nodes_[i]));
    }
}

void KnnMultiLabelTree::BuildOnRecords() {
    leaf_nodes_.resize(current_records_.size());
    for (int i = 0; i < current_records_.size(); ++i) {
        leaf_nodes_[i] = new CLeaf(mv_dataset_, current_records_[i]);
        leaf_nodes_[i]->Update(true, false);

        id_node_map_.insert(map<int, CNode*>::value_type(leaf_nodes_[i]->id(), leaf_nodes_[i]));
    }
}

