#include "tree_map_view.h"
#include <QtWidgets/QGraphicsTextItem>
#include "tree_map_item.h"
#include "variable_item.h"

TreeMapView::TreeMapView() {
	tree_item_ = NULL;
	scene_ = NULL;

	is_treemap_visible_ = false;
	is_table_lens_visible_ = false;
	highlight_var_index_ = -1;

	this->setFocusPolicy(Qt::StrongFocus);
	this->setAlignment(Qt::AlignHCenter | Qt::AlignTop);
}

TreeMapView::~TreeMapView() {

}

void TreeMapView::SetData(CNode* root, int var_num, std::vector< CNode* >& selected_nodes, int selected_count, std::vector< int >& order, std::vector< QString >& names) {
	this->root_node_ = root;
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

	if (tree_item_ == NULL) {
		tree_item_ = new TreeMapItem;
		scene_->addItem(tree_item_);
		connect(tree_item_, SIGNAL(NodeSelected(int)), this, SIGNAL(NodeSelected(int)));
	}

	if (var_items_.size() == 0) {
		var_items_.resize(var_num);
		for (int i = 0; i < var_num; ++i) {
			var_items_[i] = new VariableItem(i);
			scene_->addItem(var_items_[i]);

			connect(var_items_[i], SIGNAL(VarSelected(int)), this, SLOT(OnVarSelected(int)));
		}
	}

	var_order_ = order;

	tree_item_->SetData(root_node_);
	this->UpdateVariableItems(selected_nodes, selected_count, names);

	this->UpdateLayout();
}

void TreeMapView::SetTreeMapVisible(bool visible) {
	this->is_treemap_visible_ = visible;
	for (int i = 0; i < var_items_.size(); ++i)
		var_items_[i]->SetAbsWidthEnabled(!this->is_treemap_visible_);

	this->UpdateLayout();
}

void TreeMapView::SetTableLensVisible(bool visible) {
	this->is_table_lens_visible_ = visible;

	this->UpdateLayout();
}

void TreeMapView::UpdateLayout() {
	int total_height = 0;
	
	if (tree_item_ != NULL) {
		tree_item_->setVisible(is_treemap_visible_);
		tree_item_->setPos(0, 0);
		if (is_treemap_visible_) total_height += tree_item_->GetSize().height();
	}

	if (is_table_lens_visible_) {
		int item_margin = 10;
		total_height += item_margin;
		int var_height = 30;
		for (int i = 0; i < var_items_.size(); ++i) {
			var_items_[var_order_[i]]->setVisible(true);
			var_items_[var_order_[i]]->setPos(0, total_height);
			total_height += var_height + item_margin;
		}
	} else {
		for (int i = 0; i < var_items_.size(); ++i) {
			var_items_[var_order_[i]]->setVisible(false);
		}
	}

	QSize temp = tree_item_->GetSize();
	scene_->setSceneRect(-50, 0, temp.width(), total_height);

	scene_->update();
	this->update();
}

void TreeMapView::UpdateVariableItems(std::vector< CNode* >& selected_nodes, int selected_count, std::vector< QString >& names) {
	std::vector< std::vector< float > > var_values;
	std::vector< std::vector< int > > node_count;
	std::vector< std::vector< QColor > > node_color;
	var_values.resize(var_num_);
	node_count.resize(var_num_);
	node_color.resize(var_num_);

	for (int i = 0; i < var_num_; ++i) {
		QString var_name = names[i];
		
		for (int j = 0; j < selected_nodes.size(); ++j) {
			var_values[i].push_back(selected_nodes[j]->average_values[i]);
			node_count[i].push_back(selected_nodes[j]->point_count);
			node_color[i].push_back(selected_nodes[j]->color);
		}

		var_items_[i]->SetData(var_name, var_values[i], node_count[i], node_color[i], selected_count);
	}
}

void TreeMapView::OnVarSelected(int var_index)
{
	if (highlight_var_index_ != var_index) {
		highlight_var_index_ = var_index;
	} else {
		highlight_var_index_ = -1;
	}
	emit HighlightVarChanged(highlight_var_index_);

	for (int i = 0; i < var_items_.size(); ++i)
		var_items_[i]->SetHighlightEnabled(false);
	if (highlight_var_index_ != -1)
		var_items_[highlight_var_index_]->SetHighlightEnabled(true);
}

void TreeMapView::SetHighlightVarIndex(int index)
{
	highlight_var_index_ = index;
	for (int i = 0; i < var_items_.size(); ++i)
		var_items_[i]->SetHighlightEnabled(false);
	if (highlight_var_index_ != -1)
		var_items_[highlight_var_index_]->SetHighlightEnabled(true);
}
