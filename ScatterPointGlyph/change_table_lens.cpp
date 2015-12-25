#include "change_table_lens.h"
#include <QtGui/QPainter>
#include "change_sorting_item.h"

ChangeTableLens::ChangeTableLens() {
	scene_ = NULL;
	sorting_item_ = NULL;

	this->setFocusPolicy(Qt::StrongFocus);

	this->setAlignment(Qt::AlignLeft | Qt::AlignTop);
}

ChangeTableLens::~ChangeTableLens() {
}

void ChangeTableLens::SetData(TransMapData* data) {
	this->dataset_ = data;

	if (sorting_item_ == NULL) {
		sorting_item_ = new ChangeSortingItem;
	}
	sorting_item_->SetData(dataset_);

	if (scene_ == NULL) {
		scene_ = new QGraphicsScene(this);
		scene_->setItemIndexMethod(QGraphicsScene::NoIndex);
		setScene(scene_);

		setCacheMode(CacheBackground);
		setViewportUpdateMode(BoundingRectViewportUpdate);
		setRenderHint(QPainter::Antialiasing);
		setTransformationAnchor(AnchorUnderMouse);

		this->setBackgroundBrush(Qt::color0);
		this->autoFillBackground();

		scene_->addItem(sorting_item_);
	}

	this->update();
}

void ChangeTableLens::UpdateTableLens() {
	
}