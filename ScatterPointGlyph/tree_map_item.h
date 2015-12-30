#ifndef TREE_MAP_ITEM_H_
#define TREE_MAP_ITEM_H_

#include <QtWidgets/QGraphicsItem>
#include "tree_common.h"

#include <map>

class TreeMapItem : public QObject, public QGraphicsItem
{
	Q_OBJECT

public:
	TreeMapItem();
	~TreeMapItem();

	void SetData(CNode* data);

signals:
	void NodeSelected(int node_id);

protected:
	QRectF boundingRect() const;
	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget /* = 0 */);

	void mousePressEvent(QGraphicsSceneMouseEvent *event);

private:
	CNode* root_ = NULL;

	std::map< int, QRectF > item_pos_map_;
	std::map< int, CNode* > item_map_;

	int item_size = 15;
	int item_margin = 3;
	int transition_width = 70;
	int left_margin = 30;
	int total_width = 500;
	int total_height = 100;

	void PaintItem(QPainter* painter, CNode* node, int& bottom);
	void UpdateSize(CNode* node, int& bottom, int& item_width);
};

#endif