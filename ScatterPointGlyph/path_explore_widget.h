#ifndef PATH_EXPLORE_WIDGET_H_
#define PATH_EXPLORE_WIDGET_H_

#include <QtWidgets/QGraphicsView>
#include "path_dataset.h"

class PathExploreItem;

class PathExploreWidget : public QGraphicsView
{
	Q_OBJECT
public:
	PathExploreWidget();
	~PathExploreWidget();

	void SetData(PathDataset* data);

protected:
	void resizeEvent(QResizeEvent *event);

private:
	PathDataset* dataset_;
	
	QGraphicsScene* scene_;
	std::vector< PathExploreItem* > item_vec_;

	QGraphicsTextItem* title_item_;

	void UpdateLayout();

	private slots:
		void OnItemUpdated();
};

#endif