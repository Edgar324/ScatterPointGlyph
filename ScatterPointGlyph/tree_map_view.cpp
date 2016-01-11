#include "tree_map_view.h"
#include <QtWidgets/QGraphicsTextItem>
#include "tree_map_item.h"
#include "variable_item.h"

TreeMapView::TreeMapView() {
	tree_item_ = NULL;
	scene_ = NULL;

	this->setFocusPolicy(Qt::StrongFocus);
	this->setAlignment(Qt::AlignHCenter | Qt::AlignTop);
}

TreeMapView::~TreeMapView() {

}

void TreeMapView::SetData(CNode* data, int var_num, std::vector< CNode* >& selected_nodes, int selected_count) {
	this->root_node_ = data;
	this->var_num_ = var_num;

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

	if (var_items_.size() == 0) {
		var_items_.resize(var_num);
		for (int i = 0; i < var_num; ++i) {
			var_items_[i] = new VariableItem;
			scene_->addItem(var_items_[i]);
		}
	}

	tree_item_->SetData(root_node_);
	this->UpdateVariableItems(selected_nodes, selected_count);

	QSize temp = tree_item_->GetSize();
	int var_height = 30;
	for (int i = 0; i < var_items_.size(); ++i)
		var_items_[i]->setPos(0, (var_height + 5) * i + temp.height() + 20);

	scene_->setSceneRect(-50, 0, temp.width() + 50, temp.height() + var_items_.size() * (var_height + 5));

	scene_->update();
	this->update();
}

void TreeMapView::UpdateVariableItems(std::vector< CNode* >& selected_nodes, int selected_count) {
	std::vector< std::vector< float > > var_values;
	std::vector< std::vector< int > > node_count;
	std::vector< std::vector< QColor > > node_color;
	var_values.resize(var_num_);
	node_count.resize(var_num_);
	node_color.resize(var_num_);

	for (int i = 0; i < var_num_; ++i) {
		QString var_name = "test";
		

		for (int j = 0; j < selected_nodes.size(); ++j) {
			var_values[i].push_back(selected_nodes[j]->average_values[i]);
			node_count[i].push_back(selected_nodes[j]->point_count);
			node_color[i].push_back(selected_nodes[j]->color);

			var_items_[i]->SetData(var_name, var_values[i], node_count[i], node_color[i], selected_count);
		}
	}
}