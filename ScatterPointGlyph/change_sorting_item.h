#ifndef CHANGE_SORTING_ITEM_H_
#define CHANGE_SORTING_ITEM_H_

#include <QtWidgets/QGraphicsItem>

class TransMapData;

class ChangeSortingItem : public QObject, public QGraphicsItem
{
	Q_OBJECT

public:
	ChangeSortingItem();
	~ChangeSortingItem();

	void SetData(TransMapData* data);

protected:
	QRectF boundingRect() const;
	void mousePressEvent(QGraphicsSceneMouseEvent *event);
	void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
	void hoverMoveEvent(QGraphicsSceneHoverEvent *event);
	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget /* = 0 */);

private:
	TransMapData* trans_data_;

	std::vector< std::vector< float > > n2n_change_values_;

	int height_per_row_ = 40;
	int width_per_colume_ = 20;
	int total_height_ = 100;
	int total_width_ = 100;
	int var_label_width_ = 70;
	int item_margin_ = 10;

	void ProcessData();
};

#endif