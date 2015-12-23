#ifndef PATH_EXPLORE_ITEM_H_
#define PATH_EXPLORE_ITEM_H_

#include <QtWidgets/QGraphicsItem>
#include "path_dataset.h"

class PathExploreItem : public QObject, public QGraphicsItem
{
	Q_OBJECT

public:
	PathExploreItem();
	~PathExploreItem();

	void SetData(PathRecord* record);

	virtual QRectF boundingRect() const;
	virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

	void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event);
	void mousePressEvent(QGraphicsSceneMouseEvent *event);
	void mouseMoveEvent(QGraphicsSceneMouseEvent *event);

signals:
	void ItemUpdated();

private:
	PathRecord* path_record_;
	bool is_extending_;

	const int width_per_item_ = 50;
	const int width_per_band_ = 80;
	const int extending_height_ = 200;
	const int abstract_height_ = 50;
	const int total_width_ = 400;
	const int left_margin = 20;

	void PaintExtendingWidget(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
	void PaintAbstractWidget(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
	void PaintClusterItem(QPainter* painter, int radius, int centerx, int centery, int item_index);
	void PaintTransitionBand(QPainter* painter, int beginx, int endx, int item_index);
	QColor GetMappingColor(float value);
};

#endif