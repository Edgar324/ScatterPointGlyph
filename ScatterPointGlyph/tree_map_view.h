#ifndef TREE_MAP_VIEW_H_
#define TREE_MAP_VIEW_H_

#include <QtWidgets/QGraphicsView>
#include "tree_common.h"

class QGraphicsTextItem;
class TreeMapItem;

class TreeMapView : public QGraphicsView
{
	Q_OBJECT

public:
	TreeMapView();
	~TreeMapView();

	void SetData(CNode* data);

signals:
	void NodeSelected(int node_id);

private:
	CNode* root_node_;

	QGraphicsScene* scene_;

	TreeMapItem* tree_item_;
};

#endif