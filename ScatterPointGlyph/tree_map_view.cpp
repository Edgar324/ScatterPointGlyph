#include "tree_map_view.h"
#include "tree_map_item.h"

TreeMapView::TreeMapView() {
	tree_item_ = NULL;
	scene_ = NULL;
}

TreeMapView::~TreeMapView() {

}

void TreeMapView::SetData(CNode* data) {
	this->root_node_ = data;

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

	if (tree_item_ == NULL) {
		tree_item_ = new TreeMapItem;
		scene_->addItem(tree_item_);
	}

	tree_item_->SetData(root_node_);

	scene_->update();
	this->update();
}