/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "connector_widget.h"
#include <vtkCallbackCommand.h>

vtkStandardNewMacro(ConnectorWidget);

ConnectorWidget::ConnectorWidget() {
    this->EventCallbackCommand->SetCallback(ConnectorWidget::ProcessEvents);
}

ConnectorWidget::~ConnectorWidget() {

}

void ConnectorWidget::SetData() {

    this->BuildRepresentation();
}

void ConnectorWidget::SetEnabled(int enabled) {
    
}

void ConnectorWidget::Modified() {

}

void ConnectorWidget::BuildRepresentation() {
    
}

void ConnectorWidget::ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata) {
    ConnectorWidget* self = reinterpret_cast<ConnectorWidget*>(clientdata);

	switch (event)
	{
	case vtkCommand::LeftButtonPressEvent:
		self->OnLeftButtonDown();
		break;
	case vtkCommand::LeftButtonReleaseEvent:
		self->OnLeftButtonUp();
		break;
	case vtkCommand::RightButtonPressEvent:
		self->OnRightButtonDown();
		break;
	case vtkCommand::RightButtonReleaseEvent:
		self->OnRightButtonUp();
		break;
	case vtkCommand::MouseMoveEvent:
		self->OnMouseMove();
		break;
	default:
		break;
	}
}

void ConnectorWidget::OnMouseMove() {

}

void ConnectorWidget::OnLeftButtonDown() {

}

void ConnectorWidget::OnLeftButtonUp() {

}

void ConnectorWidget::OnRightButtonDown() {

}

void ConnectorWidget::OnRightButtonUp() {

}