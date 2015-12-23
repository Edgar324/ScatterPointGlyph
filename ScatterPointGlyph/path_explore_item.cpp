#include "path_explore_item.h"
#include <QtGui/QPainter>

PathExploreItem::PathExploreItem() {
	is_extending_ = false;
}

PathExploreItem::~PathExploreItem() {
}

void PathExploreItem::SetData(PathRecord* record) {
	this->path_record_ = record;
	this->update();
}

QRectF PathExploreItem::boundingRect() const {
	if (is_extending_) {
		return QRectF(0, 0, total_width_, extending_height_);
	} else {
		return QRectF(0, 0, total_width_, abstract_height_);
	}
}

void PathExploreItem::mousePressEvent(QGraphicsSceneMouseEvent *event) {

}

void PathExploreItem::mouseMoveEvent(QGraphicsSceneMouseEvent *event) {
}

void PathExploreItem::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) {
	is_extending_ = !is_extending_;
	
	emit ItemUpdated();
}

void PathExploreItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
	if (this->path_record_ == NULL) return;

	// paint transition band
	for (int i = 0; i < this->path_record_->change_values.size(); ++i) {
		int beginx = left_margin + i * (width_per_band_ + width_per_item_) + width_per_item_;
		int endx = beginx + width_per_band_;

		this->PaintTransitionBand(painter, beginx, endx, i);
	}

	// paint cluster item
	for (int i = 0; i < this->path_record_->item_values.size(); ++i) {
		if (!is_extending_)
			this->PaintClusterItem(painter, width_per_item_ / 2, left_margin + i * (width_per_item_ + width_per_band_) + 0.5 * width_per_item_, abstract_height_ / 2, i);
		else
			this->PaintClusterItem(painter, width_per_item_ / 2, left_margin + i * (width_per_item_ + width_per_band_) + 0.5 * width_per_item_, extending_height_ / 2, i);
	}

	// paint focus variable

	// paint highlight selection
}

void PathExploreItem::PaintAbstractWidget(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
	
}

void PathExploreItem::PaintExtendingWidget(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {

}

void PathExploreItem::PaintClusterItem(QPainter* painter, int radius, int centerx, int centery, int item_index) {
	painter->setPen(Qt::black);
	painter->drawRect(centerx - radius, centery - radius, radius * 2, radius * 2);
}

void PathExploreItem::PaintTransitionBand(QPainter* painter, int beginx, int endx, int item_index) {
	int var_num = this->path_record_->change_values[item_index].size();

	if (is_extending_) {
		for (int i = 0; i < var_num; ++i) {
			int h = extending_height_ * 0.5 - 0.5 * abstract_height_ + i * (float)abstract_height_ / var_num;
			int next_h = extending_height_ * 0.5  - 0.5 * abstract_height_ + (i + 1) * (float)abstract_height_ / var_num;
			int bandh = i * (float)extending_height_ / var_num;
			int band_nexth = (i + 1) * (float)extending_height_ / var_num;
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
			int h = i * (float)abstract_height_ / var_num;
			int next_h = (i + 1) * (float)abstract_height_ / var_num;
			painter->fillRect(QRectF(beginx, h, endx - beginx, next_h - h), this->GetMappingColor(this->path_record_->change_values[item_index][i]));
			painter->setPen(Qt::gray);
			painter->drawRect(QRectF(beginx, h, endx - beginx, next_h - h));
		}
	}
}

QColor PathExploreItem::GetMappingColor(float value) {
	if (value < 0)  return QColor(255, 0, 0);
	else return QColor(0, 255, 0);
}