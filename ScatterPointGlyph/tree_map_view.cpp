#include "tree_map_view.h"
#include <QtWidgets/QGraphicsTextItem>
#include "tree_map_item.h"

TreeMapView::TreeMapView() {
	tree_item_ = NULL;
	scene_ = NULL;

	this->setFocusPolicy(Qt::StrongFocus);
	this->setAlignment(Qt::AlignHCenter | Qt::AlignTop);
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

	/*if (title_item_ == NULL) {
		title_item_ = new QGraphicsTextItem;
		title_item_->setPlainText(QString("Tree Map"));
		scene_->addItem(title_item_);
		title_item_->setPos(this->width() / 2 - 20, 10);
	}*/

	if (tree_item_ == NULL) {
		tree_item_ = new TreeMapItem;
		scene_->addItem(tree_item_);
		tree_item_->setPos(0, 0);

		connect(tree_item_, SIGNAL(NodeSelected(int)), this, SIGNAL(NodeSelected(int)));
	}

	tree_item_->SetData(root_node_);
	QSize temp = tree_item_->GetSize();
	scene_->setSceneRect(0, 0, temp.width(), temp.height());

	scene_->update();
	this->update();
}