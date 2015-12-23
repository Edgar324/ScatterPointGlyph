#include "path_explore_widget.h"
#include "path_explore_item.h"

PathExploreWidget::PathExploreWidget() 
	: QGraphicsView() {
	scene_ = NULL;
	title_item_ = NULL;

	this->setAlignment(Qt::AlignLeft | Qt::AlignTop);
}

PathExploreWidget::~PathExploreWidget() {

}

void PathExploreWidget::SetData(PathDataset* data) {
	dataset_ = data;

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
	}

	scene_->clear();
	item_vec_.clear();

	title_item_ = new QGraphicsTextItem;
	title_item_->setPlainText(QString("Path Map"));
	scene_->addItem(title_item_);
	title_item_->setPos(this->width() / 2 - 20, 500);

	for (int i = 0; i < dataset_->path_records.size(); ++i) {
		PathExploreItem* item = new PathExploreItem();
		item->SetData(dataset_->path_records[i]);
		scene_->addItem(item);
		item_vec_.push_back(item);

		connect(item, SIGNAL(ItemUpdated()), this, SLOT(OnItemUpdated()));
	}
	this->UpdateLayout();
}

void PathExploreWidget::UpdateLayout() {
	int accu_height = 50;

	title_item_->setPos(this->width() / 2 - 20, 10);

	for (int i = 0; i < item_vec_.size(); ++i) {
		item_vec_[i]->setPos(0, accu_height);
		accu_height += item_vec_[i]->boundingRect().height() + 10;
	}

	this->update();
}

void PathExploreWidget::OnItemUpdated() {
	this->UpdateLayout();
	this->scene()->update();
}

void PathExploreWidget::resizeEvent(QResizeEvent *event) {
	QGraphicsView::resizeEvent(event);
	if (title_item_ != NULL) {
		title_item_->setPos(this->width() / 2 - 20, 10);
	}
}