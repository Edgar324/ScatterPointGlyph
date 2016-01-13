#include "variable_item.h"
#include <QtGui/QPen>
#include <QtGui/QPainter>
#include <QtWidgets/QGraphicsSceneMouseEvent>

VariableItem::VariableItem() {

}

VariableItem::~VariableItem() {

}

void VariableItem::SetData(QString var_name, std::vector< float >& var_values, std::vector< int >& node_count, std::vector< QColor >& node_color, int selected_count) {
	var_name_ = var_name;
	var_values_ = var_values;
	node_count_ = node_count;
	node_color_ = node_color;
	selected_count_ = selected_count;

	total_node_count_ = 0;
	for (int i = 0; i < node_count.size(); ++i) total_node_count_ += node_count[i];

	this->total_width = 0;
	for (int i = 0; i < node_count.size(); ++i)
		if (node_count[i] < 5)
			total_width += 0.5 * item_size + item_margin;
		else
			total_width += item_size + item_margin;

	this->prepareGeometryChange();
	this->update(QRectF(0, 0, total_width, total_height));
}

QRectF VariableItem::boundingRect() const {
	return QRectF(-150, 0, total_width + 150, total_height);
}

void VariableItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
	QPen name_pen;
	name_pen.setColor(Qt::black);
	painter->setPen(name_pen);
	painter->drawText(QRectF(-150, 0, 150, total_height), var_name_, Qt::AlignBottom | Qt::AlignHCenter);

	if (true) {
		int temp_width = 0, temp_bar_width = 0;
		for (int i = 0; i < var_values_.size(); ++i) {
			if (node_count_[i] < 5) temp_bar_width = item_size * 0.5;
			else temp_bar_width = item_size;
			if (i >= selected_count_) {
				node_color_[i].setAlpha(20);
			}
			else {
				node_color_[i].setAlpha(255);
			}
			painter->fillRect(temp_width, total_height, temp_bar_width - 1, -1 * total_height * var_values_[i], node_color_[i]);

			temp_width += temp_bar_width + item_margin;
		}
	} else {
		int temp_width = 0, temp_bar_width = 0;
		for (int i = 0; i < var_values_.size(); ++i) {
			temp_bar_width = (float)node_count_[i] / total_node_count_ * total_width;

			if (i >= selected_count_) {
				node_color_[i].setAlpha(20);
			}
			else {
				node_color_[i].setAlpha(255);
			}
			painter->fillRect(temp_width, total_height, temp_bar_width - 1, -1 * total_height * var_values_[i], node_color_[i]);

			temp_width += temp_bar_width;
		}
	}

	QPen axis_pen;
	axis_pen.setColor(Qt::gray);
	axis_pen.setWidth(2);
	painter->setPen(axis_pen);
	painter->drawLine(0, total_height, total_width, total_height);
}

