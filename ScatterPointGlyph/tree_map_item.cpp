#include "tree_map_item.h"
#include <QtGui/QPen>
#include <QtGui/QPainter>

TreeMapItem::TreeMapItem() {

}

TreeMapItem::~TreeMapItem() {

}

void TreeMapItem::SetData(CNode* data) {
	root_ = data;

	int bottom = 0;
	int item_width = 0;
	this->UpdateSize(root_, bottom, item_width);

	this->total_height = bottom;
	this->total_width = item_width + left_margin;

	this->prepareGeometryChange();
	this->update(QRectF(0, 0, total_width, total_height));
}

void TreeMapItem::UpdateSize(CNode* node, int& bottom, int& item_width) {
	if (node == NULL) return;

	int temp_bottom = bottom;
	// paint the child nodes
	if (node->type() == CNode::BRANCH && node->is_expanded) {
		CBranch* branch = dynamic_cast<CBranch*>(node);
		for (int i = 0; i < branch->linked_nodes.size(); ++i) {
			UpdateSize(branch->linked_nodes[i], bottom, item_width);
			if (i != branch->linked_nodes.size() - 1) bottom += item_margin;
		}
	}

	// paint the item glyph
	int rightx, centery;
	if (bottom - temp_bottom < item_size)  bottom += item_size;
	rightx = node->level() * (item_size + transition_width) + item_size;
	if (rightx > item_width) item_width = rightx;
}

QRectF TreeMapItem::boundingRect() const {
	return QRectF(0, 0, total_width, total_height);
}

void TreeMapItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
	item_pos_map_.clear();

	if (root_ == NULL) return;

	int bottom = 0;

	PaintItem(painter, root_, bottom);

	total_height = bottom;
}

void TreeMapItem::PaintItem(QPainter* painter, CNode* node, int& bottom) {
	if (node == NULL) return;

	int temp_bottom = bottom;
	// paint the child nodes
	std::vector< int > linked_pos;
	if (node->type() == CNode::BRANCH && node->is_expanded) {
		CBranch* branch = dynamic_cast<CBranch*>(node);
		for (int i = 0; i < branch->linked_nodes.size(); ++i) {
			float item_bottom = bottom;
			PaintItem(painter, branch->linked_nodes[i], bottom);
			linked_pos.push_back((item_bottom + bottom) / 2);
			if (i != branch->linked_nodes.size() - 1) {
				bottom += item_margin;
			}
		}
	}

	// paint the item glyph
	int leftx, centery;
	if (bottom - temp_bottom < item_size)  bottom += item_size;
	centery = (bottom + temp_bottom) / 2;
	leftx = node->level() * (item_size + transition_width) + 0.5 * item_size;

	painter->setPen(Qt::black);
	painter->drawRect(QRectF(leftx, centery - item_size / 2, item_size, item_size));

	// paint linkage
	painter->setPen(Qt::red);
	int rightx = leftx + item_size;
	int control_x1 = rightx + transition_width * 0.2;
	int control_x2 = rightx + transition_width * 0.5;
	int nextx = rightx + transition_width;
	for (int i = 0; i < linked_pos.size(); ++i) {
		QPainterPath temp_path;
		temp_path.moveTo(rightx, centery);
		temp_path.cubicTo(control_x1, centery, control_x2, linked_pos[i], nextx, linked_pos[i]);
		painter->drawPath(temp_path);
	}
}