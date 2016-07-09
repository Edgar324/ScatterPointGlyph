#include "tree_map.h"
#include <QtWidgets/QGraphicsTextItem>
#include <QtGui/QMouseEvent>
#include <QtWidgets/QToolTip>
#include "tree_map_item.h"

TreeMap::TreeMap() {
	this->setFocusPolicy(Qt::StrongFocus);
	this->setAlignment(Qt::AlignHCenter | Qt::AlignTop);
}

TreeMap::~TreeMap() {

}

void TreeMap::SetData(TreeCommon* tree) {
	this->tree_ = tree;

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
		connect(tree_item_, SIGNAL(NodeSelected(int)), this, SIGNAL(NodeSelected(int)));
	}

	tree_item_->SetData(tree_->root());
}

void TreeMap::Clear()
{
	if (scene_ != NULL) {
		scene_->clear();
		scene_->update();
	}

	tree_item_ = NULL;
	this->update();
}

void TreeMap::UpdateView() {
    if (scene_ != NULL && tree_item_ != NULL) {
        tree_item_->SetData(tree_->root());
        QSize map_size = tree_item_->GetSize();
        scene_->setSceneRect(0, 0, map_size.width(), map_size.height());
    }
    if (scene_ != NULL) this->scene_->update();
}
