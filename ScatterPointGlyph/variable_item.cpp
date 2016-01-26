#include "variable_item.h"
#include <QtGui/QPen>
#include <QtGui/QPainter>
#include <QtWidgets/QGraphicsSceneMouseEvent>

VariableItem::VariableItem(int index) {
	var_index_ = index;
	is_abs_width_on_ = false;
	is_highlight_on_ = false;
	
	ranges_[0] = 0;
	ranges_[1] = 1;
}

VariableItem::~VariableItem() {

}

void VariableItem::SetData(QString var_name, std::vector< float >& var_values, 
	std::vector< int >& node_count, std::vector< QColor >& node_color, int selected_count, 
	std::vector< std::vector< float > >& context) {

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

	
	sampled_context_data_.resize(context.size());
	for (int i = 0; i < context.size(); ++i) {
		int sample_size = 10;
		if (context[i].size() < sample_size) {
			sample_size = context[i].size();
		}
		sampled_context_data_[i].resize(sample_size);
		for (int j = 0; j < sample_size; ++j) {
			int temp = (int)((float)j / sample_size * (context[i].size() - 1));
			sampled_context_data_[i][j] = context[i][temp];
		}
	}

	ranges_[0] = 0;
	ranges_[1] = 1;

	this->prepareGeometryChange();
	this->update(QRectF(0, 0, total_width, total_height));
}

void VariableItem::SetValueRange(float min_value, float max_value)
{
	ranges_[0] = min_value;
	ranges_[1] = max_value;
}


QRectF VariableItem::boundingRect() const {
	return QRectF(-170, 0, total_width + 170, total_height);
}

void VariableItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
	QPen name_pen;
	name_pen.setColor(Qt::black);
	painter->setPen(name_pen);
	painter->drawText(QRectF(-150, 0, 150, total_height), var_name_, Qt::AlignBottom | Qt::AlignHCenter);

	if (is_highlight_on_) {
		painter->fillRect(-170, total_height / 2.0, total_height / 2.0, total_height / 2.0, Qt::red);
	}

	if (!is_abs_width_on_) {
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

			if (i >= selected_count_)
				painter->setPen(QColor(128, 128, 128, 20));
			else
				painter->setPen(QColor(128, 128, 128, 255));
			for (int j = 0; j < sampled_context_data_[i].size() - 1; ++j) {
				float x1 = temp_width + (float)(temp_bar_width - 1) * j / (sampled_context_data_[i].size() - 1);
				float y1 = total_height - total_height * sampled_context_data_[i][j];
				float x2 = temp_width + (float)(temp_bar_width - 1) * (j + 1) / (sampled_context_data_[i].size() - 1);
				float y2 = total_height - total_height * sampled_context_data_[i][j + 1];

				painter->drawLine(x1, y1, x2, y2);
			}

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

			if (i >= selected_count_)
				painter->setPen(QColor(128, 128, 128, 20));
			else
				painter->setPen(QColor(128, 128, 128, 255));
			for (int j = 0; j < sampled_context_data_[i].size() - 1; ++j) {
				float x1 = temp_width + (float)(temp_bar_width - 1) * j / (sampled_context_data_[i].size() - 1);
				float y1 = total_height - total_height * sampled_context_data_[i][j];
				float x2 = temp_width + (float)(temp_bar_width - 1) * (j + 1) / (sampled_context_data_[i].size() - 1);
				float y2 = total_height - total_height * sampled_context_data_[i][j + 1];

				painter->drawLine(x1, y1, x2, y2);
			}

			temp_width += temp_bar_width;
		}
	}

	QPen axis_pen;
	axis_pen.setColor(Qt::gray);
	axis_pen.setWidth(2);
	painter->setPen(axis_pen);
	painter->drawLine(0, total_height, total_width, total_height);
}

void VariableItem::mousePressEvent(QGraphicsSceneMouseEvent *event) {
	if (event->buttons() & Qt::RightButton)
		emit VarSelected(var_index_);
}

void VariableItem::SetAbsWidthEnabled(bool enabled)
{
	is_abs_width_on_ = enabled;

	this->update();
}

void VariableItem::SetHighlightEnabled(bool enabled)
{
	is_highlight_on_ = enabled;

	this->update();
}

QString VariableItem::GetTipString()
{
	QString tip_str = var_name_ + ": ";
	for (int i = 0; i < selected_count_; ++i)
		tip_str += QString("%0, ").arg(var_values_[i] * (ranges_[1] - ranges_[0]) + ranges_[0]);
	return tip_str;
}
