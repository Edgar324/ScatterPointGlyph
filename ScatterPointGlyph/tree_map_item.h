#ifndef TREE_MAP_ITEM_H_
#define TREE_MAP_ITEM_H_

#include <QtWidgets/QGraphicsItem>
#include "tree_common.h"

#include <map>

class TreeMapItem : public QGraphicsItem
{
public:
	TreeMapItem();
	~TreeMapItem();

	void SetData(CNode* data);

protected:
	QRectF boundingRect() const;
	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget /* = 0 */);

private:
	CNode* root_ = NULL;

	std::map< int, QRectF > item_pos_map_;

	int item_size = 20;
	int item_margin = 5;
	int transition_width = 70;
	int left_margin = 10;
	int total_width = 500;
	int total_height = 100;

	void PaintItem(QPainter* painter, CNode* node, int& bottom);
	void UpdateSize(CNode* node, int& bottom, int& item_width);
};

#endif