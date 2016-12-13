/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "table_lens.h"
#include "cnode.h"
#include "tree_common.h"
#include "multivariate_dataset.h"
#include "variable_item.h"
#include "tree_map_item.h"
#include "color_mapping_generator.h"

TableLens::TableLens() {
    this->setFocusPolicy(Qt::StrongFocus);
	this->setAlignment(Qt::AlignHCenter | Qt::AlignTop);

    scene_ = new QGraphicsScene(this);
	scene_->setItemIndexMethod(QGraphicsScene::NoIndex);
	setScene(scene_);

	setCacheMode(CacheBackground);
	setViewportUpdateMode(BoundingRectViewportUpdate);
	setRenderHint(QPainter::Antialiasing);
	setTransformationAnchor(AnchorUnderMouse);

	this->setBackgroundBrush(Qt::color0);
	this->autoFillBackground();
}

TableLens::~TableLens() {

}

void TableLens::SetData(TreeCommon* tree, vector<int>& selected_cluster_ids, vector<int>& selected_var_index) {
    cluster_tree_ = tree;
    selected_cluster_ids_ = selected_cluster_ids;
    selected_var_index_ = selected_var_index;

    this->UpdateVariableItems();
    this->UpdateLayout();
}

void TableLens::UpdateVariableItems() {
    int var_num = selected_var_index_.size();
    int cluster_num = selected_cluster_ids_.size();
    MultivariateDataset* mv_dataset = cluster_tree_->mv_dataset();

    if (var_items_.size() < var_num) {
        int current_size = var_items_.size();
		var_items_.resize(var_num);
		for (int i = current_size; i < var_num; ++i) {
			var_items_[i] = new VariableItem(i);
			scene_->addItem(var_items_[i]);
		}
	}

    vector<vector<float>> var_values;
	vector<vector<int>> node_count;
	vector<vector<float>> context_data;
	var_values.resize(var_num);
	node_count.resize(var_num);
	context_data.resize(cluster_num);

    vector<CNode*> selected_nodes;
    selected_nodes.resize(cluster_num);
    for (int i = 0; i < selected_cluster_ids_.size(); ++i)
        selected_nodes[i] = cluster_tree_->GetNode(selected_cluster_ids_[i]);

	for (int i = 0; i < var_num; ++i) {
        int var_index = selected_var_index_[i];
		QString var_name = mv_dataset->var_names()[var_index];
		
		for (int j = 0; j < cluster_num; ++j) {
			var_values[i].push_back(selected_nodes[j]->mean_values()[var_index]);
			node_count[i].push_back(selected_nodes[j]->point_count());
            // TODO: Add Context data
			/*tree_->GetNodeValues(selected_nodes[j], i, context_data[j]);
			std::sort(context_data[j].begin(), context_data[j].end());*/
		}

        /// TODO: implement this part
		/*var_items_[i]->SetData(var_name, mv_dataset->var_colors[var_index], var_values[i], node_count[i], cluster_num, context_data);
		var_items_[i]->SetValueRange(mv_dataset->original_value_ranges[var_index][0], mv_dataset->original_value_ranges[var_index][1]);*/
	}
}

void TableLens::UpdateLayout() {
    int total_height = 0;
	
	if (scene_ == NULL) return;

	int item_margin = 10;
	total_height += item_margin;
	int var_height = 30;
	for (int i = 0; i < selected_var_index_.size(); ++i) {
		var_items_[i]->setVisible(true);
		var_items_[i]->setPos(0, total_height);
		total_height += var_height + item_margin;
	}

    for (int i = selected_var_index_.size(); i < var_items_.size(); ++i) 
        var_items_[i]->setVisible(false);
    
    if (var_items_.size() > 0)
	    scene_->setSceneRect(0, 0, var_items_[0]->GetWidth(), total_height);

    this->update();
}
