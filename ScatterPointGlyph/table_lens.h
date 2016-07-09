/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef TABLE_LENS_H_
#define TABLE_LENS_H_

#include <QGraphicsView>
#include <QGraphicsScene>
#include <vector>
using namespace std;

class TreeMapItem;
class TreeCommon;
class VariableItem;

class TableLens : public QGraphicsView
{
public:
    TableLens();
    ~TableLens();

    void SetData(TreeCommon* tree, vector<int>& selected_cluster_ids, vector<int>& selected_var_index);

private:
    QGraphicsScene* scene_;

	TreeMapItem* tree_item_;
	vector<VariableItem*> var_items_;

    TreeCommon* cluster_tree_;
    vector<int> selected_cluster_ids_;
    vector<int> selected_var_index_;

	void UpdateVariableItems();
	void UpdateLayout();
};

#endif