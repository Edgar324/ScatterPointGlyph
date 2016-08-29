/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "cluster_projection_tree.h"
#include "scatter_point_dataset.h"
#include "multi_label_processor.h"
#include "tsne_projector.h"
#include "mds_projector.h"

ClusterProjectionTree::ClusterProjectionTree(ScatterPointDataset* data)
    : MultiLabelTree(data) {
    leaf_nodes = root_->linked_nodes;
}

ClusterProjectionTree::~ClusterProjectionTree() {

}

void ClusterProjectionTree::ConstructTree(float left, float right, float bottom, float top, float glyph_radius) {
    // Step one: reconstruct the root node and separate leaf nodes into two branch nodes.
    // One contains leaf nodes that are in the viewport, the other contains other leaf nodes.
    for (int i = 0; i < root_->linked_nodes.size(); ++i) {
        if (root_->linked_nodes[i]->type() == CNode::BRANCH) {
            CBranch* branch = (CBranch*)(root_->linked_nodes[i]);
            branch->linked_nodes.clear();
            id_node_map_.erase(branch->id());
            delete branch;
        }
    }
    root_->linked_nodes = leaf_nodes;

    // Step two: clustering the node
    vector<vector<float>> pos;
	vector<vector<float>> value;
	for (int i = 0; i < root_->linked_nodes.size(); ++i) {
		pos.push_back(root_->linked_nodes[i]->mean_pos);
		value.push_back(root_->linked_nodes[i]->mean_values);
	}

    int point_num = root_->linked_nodes.size();
    vector<vector<bool>> connecting_status;
    connecting_status.resize(point_num);
    for (int i = 0; i < point_num; ++i)
        connecting_status[i].resize(point_num, true);

    bool is_splitted = false;
    float temp_radius = 0.15;

	// as long as the node is not split, split again with a smaller radius for the label estimation
    vector<int> clusters;
    while (!is_splitted) {
        processor_->SetLabelEstimationRadius(temp_radius);
	    processor_->SetData(pos, value, point_dataset_->var_weights, connecting_status);
	    //this->GenerateSegmentUncertainty(node->linked_nodes, connecting_status, edge_weights);
	    //processor_->SetEdgeWeights(edge_weights);
	    processor_->GenerateCluster();
        clusters = processor_->GetResultLabel();

        int first_label = -1, second_label = -1;
        for (int i = 0; i < clusters.size(); ++i) 
            if (clusters[i] >= 0) {
                first_label = clusters[i];
                for (int j = i + 1; j < clusters.size(); ++j) 
                    if (clusters[j] >= 0 && clusters[j] != first_label) {
                        second_label = clusters[j];
                        break;
                    }
                break;
            }
        if (first_label != -1 && second_label != -1 && first_label != second_label) {
            is_splitted = true;
        }
        else {
            temp_radius /= factor_;
        }
    }

    // add cluster nodes
    vector<CBranch*> label_nodes;
	label_nodes.resize(root_->linked_nodes.size(), NULL);
	for (int i = 0; i < clusters.size(); ++i) {
		int label = clusters[i];
		if (label < 0) continue;
		if (label_nodes[label] == NULL) {
			label_nodes[label] = new CBranch;
			label_nodes[label]->set_level(root_->level() + 1);
			label_nodes[label]->radius = temp_radius;
			label_nodes[label]->parent = root_;
		}
		label_nodes[label]->linked_nodes.push_back(root_->linked_nodes[i]);
		root_->linked_nodes[i]->set_level(root_->level() + 2);
		root_->linked_nodes[i]->parent = label_nodes[label];
	}

	vector<CNode*> linked_nodes;
	for (int i = 0; i < clusters.size(); ++i)
		if (label_nodes[i] != NULL) {
			linked_nodes.push_back(label_nodes[i]);
		}
	for (int i = 0; i < linked_nodes.size(); ++i)
		ProgressNodeData(linked_nodes[i]);

	root_->linked_nodes.clear();
    for (int i = 0; i < linked_nodes.size(); ++i) {
        root_->linked_nodes.push_back(linked_nodes[i]);
    }

    // project the nodes
    vector<vector<float>> center_values;
    for (int i = 0; i < root_->linked_nodes.size(); i++)
        center_values.push_back(root_->linked_nodes[i]->mean_values);
    vector<vector<float>> projected_pos;
    TsneProjector projector;
    projector.Project(center_values, projected_pos);
    for (int i = 0; i < root_->linked_nodes.size(); i++)
        root_->linked_nodes[i]->mean_pos = vector<float>{projected_pos[0][i], projected_pos[1][i]};
}
