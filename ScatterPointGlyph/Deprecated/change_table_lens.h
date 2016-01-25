#ifndef CHANGE_TABLE_LENS_H_
#define CHANGE_TABLE_LENS_H_

#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QWidget>

class TransMapData;
class ChangeSortingItem;

class ChangeTableLens : public QGraphicsView
{
public:
	ChangeTableLens();
	~ChangeTableLens();

	void SetData(TransMapData* data);

private:
	TransMapData* dataset_;

	QGraphicsScene* scene_;
	ChangeSortingItem* sorting_item_;

	void UpdateTableLens();
};

#endif