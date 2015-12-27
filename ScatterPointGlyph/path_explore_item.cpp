#include "path_explore_item.h"
#include <QtGui/QPainter>
#include <QtWidgets/QGraphicsSceneMouseEvent>
#include <QtWidgets/QToolTip>

PathExploreItem::PathExploreItem() {
	is_extending_ = false;

	this->setAcceptHoverEvents(true);

	setFlags(ItemIsSelectable);
}

PathExploreItem::~PathExploreItem() {
}

void PathExploreItem::SetData(PathRecord* record) {
	this->path_record_ = record;

	this->total_height = (size_per_item_ + row_margin_) * path_record_->change_values.size() / (item_num_per_row_ - 1);

	this->update();
}

void PathExploreItem::SetItemWidth(int w) {
	this->total_width_ = w;
	
	if (this->path_record_ != NULL) {
		this->total_height = (size_per_item_ + row_margin_) * path_record_->change_values.size() / (item_num_per_row_ - 1);
	}

	this->update();
}

void PathExploreItem::SetSelectedVar(int index) {
	this->selected_var_ = index;

	this->update();
}

QRectF PathExploreItem::boundingRect() const {
	return QRectF(0, -5, total_width_, total_height + 5);
}

void PathExploreItem::mousePressEvent(QGraphicsSceneMouseEvent *event) {
	if (this->path_record_ == NULL) return;

	QGraphicsItem::mousePressEvent(event);
}

void PathExploreItem::hoverMoveEvent(QGraphicsSceneHoverEvent *event) {
	if (this->path_record_ == NULL) return;

	int x = event->pos().x();
	int y = event->pos().y();

	int row = y / (size_per_item_ + row_margin_);
	int colume = (x - label_width_ - 2 * item_margin_) / (size_per_item_ + width_per_band_);
	int x_colume = x - label_width_ - 2 * item_margin_ - colume * (size_per_item_ + width_per_band_) - size_per_item_;
	if (x_colume > 7 && x_colume < width_per_band_ - 12) {
		float width_per_var = (float)(width_per_band_ - 19) / this->path_record_->change_values[0].size();
		int temp_index = (x_colume - 7) / width_per_var;
		if (temp_index != selected_var_) {
			selected_var_ = temp_index;
			emit SelectedVarChanged(selected_var_);
		}
	} else {
		emit SelectedVarChanged(-1);
	}

	QGraphicsItem::hoverMoveEvent(event);
}

void PathExploreItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
	if (this->path_record_ == NULL) return;

	width_per_band_ = (this->total_width_ - 4 * item_margin_ - 2 * label_width_ - item_num_per_row_ * size_per_item_) / (item_num_per_row_ - 1);
	// paint begin label
	QPen path_pen;
	path_pen.setColor(Qt::gray);
	path_pen.setWidth(2.0);

	painter->setPen(path_pen);
	painter->drawEllipse(item_margin_, size_per_item_ / 2 - label_width_ / 2, label_width_, label_width_);
	painter->drawLine(item_margin_ + label_width_, size_per_item_ / 2, item_margin_ * 2 + label_width_, size_per_item_ / 2);

	int left_item_begin = item_margin_ * 2 + label_width_;
	// paint transition band
	for (int i = 0; i < this->path_record_->change_values.size(); ++i) {
		int colume = i % (item_num_per_row_ - 1);
		int beginx = left_item_begin + colume * (width_per_band_ + size_per_item_) + size_per_item_;
		int endx = beginx + width_per_band_;

		this->PaintTransitionBand(painter, beginx, endx, i);
	}

	// paint cluster item
	for (int i = 0; i < this->path_record_->item_values.size(); ++i) {
		int row = i / (item_num_per_row_ - 1);
		int colume = i % (item_num_per_row_ - 1);
		int centerx, centery;
		if (colume == 0 && row != 0) {
			centerx = left_item_begin + (item_num_per_row_ - 1) * (size_per_item_ + width_per_band_) + 0.5 * size_per_item_;
			centery = (row - 1) * (size_per_item_ + row_margin_) + size_per_item_ / 2;
			this->PaintClusterItem(painter, size_per_item_ / 2, centerx, centery, i - 1);
			
			if (i != path_record_->item_values.size() - 1) {
				int nextx = left_item_begin + size_per_item_ / 2;
				int nexty = row * (size_per_item_ + row_margin_);
				QPainterPath path;
				path.moveTo(centerx, nexty - row_margin_);
				path.cubicTo(centerx, nexty - row_margin_ / 2, centerx, nexty - row_margin_ / 2, centerx - size_per_item_ / 2, nexty - row_margin_ / 2);
				path.lineTo(nextx + row_margin_ / 2, nexty - row_margin_ / 2);
				path.cubicTo(nextx, nexty - row_margin_ / 2, nextx, nexty - row_margin_ / 2, nextx, nexty);
				painter->setPen(path_pen);
				painter->drawPath(path);
			}
		}

		if (i != path_record_->item_values.size() - 1 || (i == path_record_->item_values.size() - 1 && colume != 0)) {
			centerx = left_item_begin + colume * (size_per_item_ + width_per_band_) + 0.5 * size_per_item_;
			centery = row * (size_per_item_ + row_margin_) + size_per_item_ / 2;
			this->PaintClusterItem(painter, size_per_item_ / 2, centerx, centery, i);
		}

		if (i == this->path_record_->item_values.size() - 1) {
			painter->setPen(path_pen);
			painter->drawEllipse(centerx + item_margin_ + size_per_item_ / 2, centery - label_width_ / 2, label_width_, label_width_);
			painter->drawLine(centerx + size_per_item_ / 2, centery, centerx + size_per_item_ / 2 + item_margin_, centery);
		}
	}

	// paint focus variable

	// paint highlight selection
}

void PathExploreItem::PaintClusterItem(QPainter* painter, int radius, int centerx, int centery, int item_index) {
	painter->setPen(Qt::gray);
	painter->drawRect(centerx - radius, centery - radius, radius * 2, radius * 2);

	// draw axis
	QPen axis_pen;
	axis_pen.setColor(Qt::gray);
	axis_pen.setWidth(2.0);
	painter->setPen(axis_pen);

	std::vector< float > x_vec;
	std::vector< float > y_vec;

	for (int i = 0; i < path_record_->item_values[item_index].size(); ++i) {
		float end_arc = i * 3.14159 * 2 / path_record_->item_values[item_index].size();

		float temp_radius = radius * 0.8;
		float x = temp_radius * cos(end_arc);
		float y = temp_radius * sin(end_arc);

		x_vec.push_back(x * path_record_->item_values[item_index][i]);
		y_vec.push_back(y * path_record_->item_values[item_index][i]);

		painter->drawLine(centerx, centery, centerx + x, centery + y);
	}

	axis_pen.setColor(path_record_->item_color[item_index]);
	painter->setPen(axis_pen);
	for (int i = 0; i < path_record_->item_values[item_index].size(); ++i) {
		painter->drawLine(centerx + x_vec[i], centery + y_vec[i], centerx + x_vec[(i + 1) % path_record_->item_values[item_index].size()], centery + y_vec[(i + 1) % path_record_->item_values[item_index].size()]);
	}
}

void PathExploreItem::PaintTransitionBand(QPainter* painter, int beginx, int endx, int item_index) {
	int var_num = this->path_record_->change_values[item_index].size();

	int centery = item_index / (item_num_per_row_ - 1) * (size_per_item_ + row_margin_) + size_per_item_ / 2;

	QPen arrow_pen;
	arrow_pen.setColor(Qt::gray);
	arrow_pen.setWidth(2.0);
	painter->setPen(arrow_pen);
	painter->drawLine(beginx + 2, centery, endx - 2, centery);
	painter->drawLine(endx - 2 - 5, centery - 5, endx - 2, centery);
	painter->drawLine(endx - 2 - 5, centery + 5, endx - 2, centery);

	QPen path_pen;
	path_pen.setColor(Qt::gray);
	path_pen.setWidth(1.0);

	QPen highlight_pen;
	highlight_pen.setColor(Qt::gray);
	highlight_pen.setWidth(1.0);
	highlight_pen.setDashOffset(2);

	float width_per_var = (float)(endx - beginx - 19) / var_num;
	for (int i = 0; i < var_num; ++i) {
		float value = path_record_->change_values[item_index][i];
		int tempx = beginx + 5 +  width_per_var * i + 2;
		int temp_height = value * size_per_item_ / 2 * -1 * 0.8;
		if (temp_height > 0){
			if (selected_var_ == -1 || i == selected_var_)
				painter->fillRect(QRectF(tempx, centery, width_per_var, temp_height), QColor(255, 0, 0, 255));
			else
				painter->fillRect(QRectF(tempx, centery, width_per_var, temp_height), QColor(255, 0, 0, 50));
		}
		else {
			if (selected_var_ == -1 || i == selected_var_)
				painter->fillRect(QRectF(tempx, centery, width_per_var, temp_height), QColor(0, 255, 0, 255));
			else
				painter->fillRect(QRectF(tempx, centery, width_per_var, temp_height), QColor(0, 255, 0, 50));
		}

		/*painter->setPen(path_pen);
		painter->drawRect(QRectF(tempx, centery, width_per_var, temp_height));*/

		if (i == selected_var_) {
			painter->setPen(highlight_pen);
			painter->drawLine(tempx, centery - size_per_item_ / 2 * 0.8, tempx + width_per_var, centery - size_per_item_ / 2 * 0.8);
			painter->drawLine(tempx, centery + size_per_item_ / 2 * 0.8, tempx + width_per_var, centery + size_per_item_ / 2 * 0.8);
			painter->drawLine(tempx, centery - size_per_item_ / 2 * 0.8, tempx, centery + size_per_item_ / 2 * 0.8);
			painter->drawLine(tempx + width_per_var, centery - size_per_item_ / 2 * 0.8, tempx + width_per_var, centery + size_per_item_ / 2 * 0.8);

			painter->setPen(Qt::black);
			QString str = QString::fromLocal8Bit(this->path_record_->var_names[selected_var_].c_str());
			painter->drawText(QPoint(beginx + 7, centery - size_per_item_ * 0.4), str + QString(": %0").arg(this->path_record_->change_values[item_index][selected_var_], 5, 'g', 3));
		}
	}
	/*if (is_extending_) {
		for (int i = 0; i < var_num; ++i) {
			int h = height_per_item_ * 0.5 - 0.5 * height_per_item_ + i * (float)height_per_item_ / var_num;
			int next_h = height_per_item_ * 0.5  - 0.5 * height_per_item_ + (i + 1) * (float)height_per_item_ / var_num;
			int bandh = i * (float)height_per_item_ / var_num;
			int band_nexth = (i + 1) * (float)height_per_item_ / var_num;
			int controlx1 = beginx + (endx - beginx) * 0.1;
			int controlx2 = beginx + (endx - beginx) * 0.9;
			int curvex1 = beginx + (endx - beginx) * 0.2;
			int curvex2 = beginx + (endx - beginx) * 0.8;
			QPainterPath path;
			path.moveTo(beginx, h);
			path.cubicTo(controlx1, h, controlx1, bandh, curvex1, bandh);
			path.lineTo(curvex2, bandh);
			path.cubicTo(controlx2, bandh, controlx2, h, endx, h);
			path.lineTo(endx, next_h);
			path.cubicTo(controlx2, next_h, controlx2, band_nexth, curvex2, band_nexth);
			path.lineTo(curvex1, band_nexth);
			path.cubicTo(controlx1, band_nexth, controlx1, next_h, beginx, next_h);
			
			painter->fillPath(path, this->GetMappingColor(this->path_record_->change_values[item_index][i]));

			painter->setPen(Qt::gray);
			painter->drawPath(path);
		}
	} else {
		for (int i = 0; i < var_num; ++i) {
			int h = i * (float)height_per_item_ / var_num;
			int next_h = (i + 1) * (float)height_per_item_ / var_num;
			painter->fillRect(QRectF(beginx, h, endx - beginx, next_h - h), this->GetMappingColor(this->path_record_->change_values[item_index][i]));
			painter->setPen(Qt::gray);
			painter->drawRect(QRectF(beginx, h, endx - beginx, next_h - h));
		}
	}*/
}

QColor PathExploreItem::GetMappingColor(float value) {
	if (value < 0)  return QColor(255, 0, 0);
	else return QColor(0, 255, 0);
}