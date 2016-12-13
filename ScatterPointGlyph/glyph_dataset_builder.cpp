/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "glyph_dataset_builder.h"
#include "multivariate_dataset.h"
#include "tree_common.h"
#include "glyph_object.h"
#include "glyph_dataset.h"
#include "utility.h"

void GlyphDatasetBuilder::Build(MultivariateDataset* mv_dataset, vector<int>& selected_var_index,
    TreeCommon* tree, float left, float right, float bottom, float top, 
    float glyph_half_size, GlyphDataset* glyph_dataset) {

    // Get visible nodes
    vector<CNode*> visible_nodes;
    tree->GetNodes(left, right, bottom, top, visible_nodes);

    // Generate axis order
    vector<vector<float>> mean_values;
    for (int i = 0; i < visible_nodes.size(); ++i) {
        vector<float> temp_mean;
        for (int j = 0; j < selected_var_index.size(); ++j)
            temp_mean.push_back(visible_nodes[i]->mean_values()[selected_var_index[j]]);
        mean_values.push_back(temp_mean);
    }
        
    vector<int> axis_order;
    axis_order.resize(selected_var_index.size());
    for (int i = 0; i < selected_var_index.size(); ++i) axis_order[i] = i;

    if (mean_values.size() > 2) {
        Utility::GenerateAxisOrder(mean_values, axis_order);
    }
    vector<int> selected_vars;
    selected_vars.resize(axis_order.size());
    for (int i = 0; i < axis_order.size(); ++i)
        selected_vars[i] = selected_var_index[axis_order[i]];

    selected_var_index = selected_vars;

    vector<QString> names;
    vector<QColor> colors;
    for (int i = 0; i < selected_vars.size(); ++i) {
        names.push_back(mv_dataset->var_names()[selected_vars[i]]);
        colors.push_back(mv_dataset->var_colors()[selected_vars[i]]);
    }

    int max_point_count = -1;
    for (int i = 0; i < visible_nodes.size(); ++i)
        if (visible_nodes[i]->point_count() > max_point_count)
            max_point_count = visible_nodes[i]->point_count();

    // NOTE: Clear and call Modified() before updating the rendering view
    vector<float> node_saliency;
    EvaluateSaliency(mv_dataset, visible_nodes, node_saliency);

    glyph_dataset->Clear();
    for (int i = 0; i < visible_nodes.size(); ++i) {
        CNode* node = visible_nodes[i];

        vector<float> means;
        vector<float> std_devs;

        for (int j = 0; j < selected_vars.size(); ++j) {
            int var_index = selected_vars[j];
            means.push_back(node->mean_values()[var_index]);
            std_devs.push_back(node->std_deviations()[var_index]);
        }

        bool is_expandable = false;
        if (node->type() == CNode::BRANCH) {
            is_expandable = true;
            CBranch* branch = (CBranch*)node;
            for (int j = 0; j < branch->children().size(); ++j) {
                if (branch->children()[j]->type() == CNode::LEAF) {
                    is_expandable = false;
                    break;
                }
            }
        }

        GlyphObject* object = new GlyphObject(node->id(),
            names, colors, means, std_devs, node_saliency[i], node->point_count(), max_point_count, 
            glyph_half_size, node->mean_pos()[0], node->mean_pos()[1], is_expandable);

        glyph_dataset->AddGlyphObject(object);
    }

    // TODO: construct the border lines for each node
}

void GlyphDatasetBuilder::EvaluateSaliency(MultivariateDataset* mv_dataset, vector<CNode*>& nodes, vector<float>& node_saliency) {
    node_saliency.clear();
    node_saliency.resize(nodes.size(), 0.5);

    //if (nodes.size() > 1) {
    //    vector<vector<bool>> is_connecting;
    //    float min_edge_length;
    //    Utility::Triangulation(nodes, is_connecting, min_edge_length);

    //    // update node saliency
    //    node_saliency.resize(nodes.size());
    //    node_saliency.assign(nodes.size(), 0.0);
    //    for (int i = 0; i < nodes.size(); ++i) {
    //        for (int j = 0; j < nodes.size(); ++j)
    //            if (i != j && is_connecting[i][j]) {
    //                float var_dis = 0;
    //                for (int k = 0; k < mv_dataset->var_num(); ++k)
    //                    var_dis += abs(nodes[i]->mean_values()[k] - nodes[j]->mean_values()[k]) * mv_dataset->var_weights()[k];
    //                float space_dis = 0;
    //                space_dis = sqrt(pow(nodes[i]->mean_pos()[0] - nodes[j]->mean_pos()[0], 2) + pow(nodes[i]->mean_pos()[1] - nodes[j]->mean_pos()[1], 2));
    //                node_saliency[i] += var_dis * exp(-1 * space_dis / 0.25);
    //            }
    //    }
    //    float max_saliency = -1e10;
    //    for (int i = 0; i < nodes.size(); ++i)
    //        if (max_saliency < node_saliency[i]) max_saliency = node_saliency[i];

    //    if (max_saliency > 0) {
    //        for (int i = 0; i < nodes.size(); ++i)
    //            node_saliency[i] /= max_saliency;
    //    } else {
    //        for (int i = 0; i < nodes.size(); ++i)
    //            node_saliency[i] = 0.5;
    //    }
    //}
    //else {
    //    node_saliency.resize(nodes.size());
    //    node_saliency.assign(nodes.size(), 0.5);
    //}
}
