#ifndef TREE_MAP_H_
#define TREE_MAP_H_

#include <QtWidgets/QGraphicsView>
#include "tree_common.h"

class QGraphicsTextItem;
class TreeMapItem;
class VariableItem;

class TreeMap : public QGraphicsView
{
	Q_OBJECT

public:
	TreeMap();
	~TreeMap();

	void SetData(TreeCommon* tree);
    void UpdateView();
	void Clear();

signals:
	void NodeSelected(int node_id);

private:
	TreeCommon* tree_ = NULL;
	TreeMapItem* tree_item_ = NULL;

    QGraphicsScene* scene_ = NULL;
};

#endif