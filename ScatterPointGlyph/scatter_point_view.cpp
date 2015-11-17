#include "scatter_point_view.h"
#include <QtCore/QTimer>
#include <QtGui/QWheelEvent>
#include <vtkRenderWindow.h>

ScatterPointView::ScatterPointView() 
	: is_wheel_updated_(false) {
	timer_ = new QTimer;
	timer_->setSingleShot(true);

	connect(timer_, SIGNAL(timeout()), this, SLOT(OnTimeout()));
}

ScatterPointView::~ScatterPointView() {

}

void ScatterPointView::wheelEvent(QWheelEvent* event) {
	is_wheel_updated_ = true;
	if (!timer_->isActive()) {
		timer_->start(500);
	}

	vtkRenderWindowInteractor* interactor = NULL;

	if (this->mRenWin) {
		interactor = this->mRenWin->GetInteractor();
	}

	if (!interactor || !interactor->GetEnabled()) {
		return;
	}

	if (event->delta() > 0)
		interactor->InvokeEvent(vtkCommand::MouseWheelForwardEvent, event);
	else
		interactor->InvokeEvent(vtkCommand::MouseWheelBackwardEvent, event);
}

void ScatterPointView::OnTimeout() {
	if (is_wheel_updated_) {
		is_wheel_updated_ = false;
		timer_->start(500);
	} else {
		emit ViewUpdated();
	}
}