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
	void SetItemWidth(int w);
	void SetSelectedVar(int index);

	virtual QRectF boundingRect() const;
	virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

	void mousePressEvent(QGraphicsSceneMouseEvent *event);
	void hoverMoveEvent(QGraphicsSceneHoverEvent *event);
	//void mouseMoveEvent(QGraphicsSceneMouseEvent *event);

signals:
	void ItemUpdated();
	void SelectedVarChanged(int);

private:
	PathRecord* path_record_;
	bool is_extending_;

	int item_num_per_row_ = 3;
	int width_per_item_ = 50;
	int width_per_band_ = 80;
	int height_per_item_ = 50;
	int total_width_ = 400;
	int total_height = 100;
	int row_margin_ = 15;

	int item_margin_ = 10;
	int label_width_ = 10;

	int selected_var_ = -1;


	void PaintClusterItem(QPainter* painter, int radius, int centerx, int centery, int item_index);
	void PaintTransitionBand(QPainter* painter, int beginx, int endx, int item_index);
	QColor GetMappingColor(float value);
};

#endif