#ifndef TREE_MAP_VIEW_H_
#define TREE_MAP_VIEW_H_

#include <QtWidgets/QGraphicsView>
#include "tree_common.h"

class QGraphicsTextItem;
class TreeMapItem;
class VariableItem;

class TreeMapView : public QGraphicsView
{
	Q_OBJECT

public:
	TreeMapView();
	~TreeMapView();

	void SetData(CNode* data, int var_num, std::vector< CNode* >& selected_nodes, int selected_count, std::vector< int >& order, std::vector< QString >& names);

signals:
	void NodeSelected(int node_id);

private:
	CNode* root_node_;
	int var_num_;
	std::vector< int > var_order_;

	QGraphicsScene* scene_;

	TreeMapItem* tree_item_;
	std::vector< VariableItem* > var_items_;

	void UpdateVariableItems(std::vector< CNode* >& selected_nodes, int selected_count, std::vector< QString >& names);
};

#endif