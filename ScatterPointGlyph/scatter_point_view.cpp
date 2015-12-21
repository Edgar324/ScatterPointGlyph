#include "scatter_point_view.h"
#include <iostream>
#include <QtCore/QTimer>
#include <QtGui/QWheelEvent>
#include <QtGui/QDragMoveEvent>
#include <vtkRenderWindow.h>

ScatterPointView::ScatterPointView() 
	: QVTKWidget(), is_wheel_updated_(false) {
	timer_ = new QTimer;
	timer_->setSingleShot(true);

	connect(timer_, SIGNAL(timeout()), this, SLOT(OnTimeout()));
}

ScatterPointView::~ScatterPointView() {

}

void ScatterPointView::wheelEvent(QWheelEvent* event) {
	is_wheel_updated_ = true;
	if (!timer_->isActive()) {
		timer_->start(200);
	}

	QVTKWidget::wheelEvent(event);
}

void ScatterPointView::mousePressEvent(QMouseEvent* event) {
	QVTKWidget::mousePressEvent(event);
	emit GlyphSelected(event->x(), event->y());
}

void ScatterPointView::mouseMoveEvent(QMouseEvent* event) {
	if (event->buttons() & Qt::MidButton) QVTKWidget::mouseMoveEvent(event);
}

void ScatterPointView::OnTimeout() {
	if (is_wheel_updated_) {
		is_wheel_updated_ = false;
		timer_->start(200);
	} else {
		emit ViewUpdated();
	}
}