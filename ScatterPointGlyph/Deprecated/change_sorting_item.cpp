#include "change_sorting_item.h"
#include "transmap_data.h"
#include <QtGui/QPen>
#include <QtGui/QPainter>

ChangeSortingItem::ChangeSortingItem() {
	this->setAcceptHoverEvents(true);
}

ChangeSortingItem::~ChangeSortingItem() {

}

void ChangeSortingItem::SetData(TransMapData* data) {
	trans_data_ = data;

	int node_num = data->level_one_nodes.size();
	this->total_width_ = 2 * item_margin_ + var_label_width_ + width_per_colume_ * node_num * (node_num - 1);
	this->total_height_ = data->var_num * height_per_row_;

	ProcessData();
}

void ChangeSortingItem::ProcessData() {
	int node_num = trans_data_->level_one_nodes.size();
	n2n_change_values_.clear();
	for (int i = 0; i < node_num; ++i)
		for (int j = 0; j < node_num; ++j)
			if (j != i) {
				std::vector< float > temp_value;
				temp_value.resize(trans_data_->var_num);
				for (int k = 0; k < trans_data_->var_num; ++k)
					temp_value[k] = trans_data_->level_one_nodes[i]->average_values[k] - trans_data_->level_one_nodes[j]->average_values[k];
				n2n_change_values_.push_back(temp_value);
			}
}

QRectF ChangeSortingItem::boundingRect() const {
	return QRectF(0, 0, this->total_width_, this->total_height_);
}

void ChangeSortingItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget ) {
	// paint variable label
	QPen label_pen;
	label_pen.setColor(Qt::black);
	painter->setPen(label_pen);

	for (int i = 0; i < trans_data_->var_num; ++i) {
		int w = item_margin_;
		int h = i * height_per_row_;
		
		painter->drawText(w, h, var_label_width_, height_per_row_, Qt::AlignLeft | Qt::AlignVCenter, "V");
	}

	QPen axis_pen;
	axis_pen.setColor(Qt::gray);
	axis_pen.setWidth(2.0);

	// paint variable center line and variable change values
	int normalized_height = height_per_row_ / 2 - 2;

	for (int i = 0; i < trans_data_->var_num; ++i) {
		int centerh = i * height_per_row_ + 0.5 * height_per_row_;

		for (int j = 0; j < n2n_change_values_.size(); ++j) {
			int beginx = 2 * item_margin_ + var_label_width_ + width_per_colume_ * j;
			int endx = beginx + width_per_colume_ ;

			if (n2n_change_values_[j][i] > 0) {
				painter->fillRect(beginx, centerh, width_per_colume_, -1 * normalized_height * n2n_change_values_[j][i], Qt::green);
			} else {
				painter->fillRect(beginx, centerh, width_per_colume_, -1 * normalized_height * n2n_change_values_[j][i], Qt::red);
			}
			painter->drawRect(beginx, centerh, width_per_colume_, -1 * normalized_height * n2n_change_values_[j][i]);
		}
	}
}

void ChangeSortingItem::mousePressEvent(QGraphicsSceneMouseEvent *event) {

}

void ChangeSortingItem::hoverMoveEvent(QGraphicsSceneHoverEvent *event) {

}

void ChangeSortingItem::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {

}
