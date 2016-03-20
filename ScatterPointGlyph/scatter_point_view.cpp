#include "scatter_point_view.h"
#include <iostream>
#include <QtCore/QTimer>
#include <QtGui/QWheelEvent>
#include <QtGui/QDragMoveEvent>
#include <QtWidgets/QToolTip>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleImage.h>
#include "cnode.h"

ScatterPointView::ScatterPointView() 
	: QVTKWidget(), is_wheel_updated_(false) {
	timer_ = new QTimer;
	timer_->setSingleShot(true);

	this->setFocusPolicy(Qt::StrongFocus);

	vtkSmartPointer< vtkInteractorStyleImage > imageStyle =
		vtkSmartPointer< vtkInteractorStyleImage >::New();
	this->GetRenderWindow()->GetInteractor()->SetInteractorStyle(imageStyle);
	this->GetRenderWindow()->PointSmoothingOff();

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
	if (event->buttons() & Qt::LeftButton)
		emit GlyphSelected(event->x(), this->height() - event->y());
	if (event->buttons() & Qt::RightButton)
		emit RightButtonDown();
}

void ScatterPointView::mouseMoveEvent(QMouseEvent* event) {
	global_mouse_pos_ = event->globalPos();
	if (event->buttons() & Qt::LeftButton) {
		emit MouseDrag(event->pos().x(), this->height() - event->pos().y());
	}
	QVTKWidget::mouseMoveEvent(event);
}

void ScatterPointView::mouseReleaseEvent(QMouseEvent* event) {
	QVTKWidget::mouseReleaseEvent(event);
	if (event->button() == Qt::LeftButton)
		emit LeftButtonUp();
}

void ScatterPointView::OnTimeout() {
	if (is_wheel_updated_) {
		is_wheel_updated_ = false;
		timer_->start(200);
	} else {
		emit ViewUpdated();
	}
}

void ScatterPointView::ShowTooltip(int point_count, QString axis_name, float average_value, float variance_value)
{
	QString tip_str = "";
	tip_str = "Cluster Information: <br>";
	tip_str += QString("Point Count: %0<br>").arg(point_count);
	tip_str += "Variable Values: Mean(SD)<br>";
	QString temp_str = "";
	temp_str += "<font color=red>";
	temp_str += axis_name;
	temp_str += QString(":  %0(%1)").arg(average_value).arg(variance_value);
	temp_str += "</font><br>";
	tip_str += temp_str;

	QToolTip::showText(global_mouse_pos_, tip_str);
}

void ScatterPointView::HideTooltip()
{
	QToolTip::hideText();
}

void ScatterPointView::SetHighlightVarIndex(int var_index)
{
    if (var_index == -1) QToolTip::hideText();
	emit HighlightVarChanged(var_index);
}
