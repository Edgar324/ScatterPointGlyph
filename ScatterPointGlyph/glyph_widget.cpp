/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "glyph_widget.h"
#include "glyph_object.h"
#include <vtkCallbackCommand.h>

vtkStandardNewMacro(GlyphWidget);

GlyphWidget::GlyphWidget() {
    this->EventCallbackCommand->SetCallback(GlyphWidget::ProcessEvents);
}

GlyphWidget::~GlyphWidget() {

}

void GlyphWidget::SetData(GlyphObject* object) {
    this->glyph_object_ = object;
    if (this->glyph_object_ != NULL) 
        this->glyph_object_->set_glyph_widget(this);
    this->BuildRepresentation();
}

void GlyphWidget::SetRenderingWidget(GlyphRenderingWidget* widget) {
    this->rendering_widget_ = widget;
}

void GlyphWidget::SetEnabled(int enabled) {
    
}

void GlyphWidget::BuildRepresentation() {
    
}

void GlyphWidget::ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata) {
    GlyphWidget* self = reinterpret_cast<GlyphWidget*>(clientdata);

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

void GlyphWidget::OnMouseMove() {
    cout << "Glyph widget move" << endl;
}

void GlyphWidget::OnLeftButtonDown() {

}

void GlyphWidget::OnLeftButtonUp() {

}

void GlyphWidget::OnRightButtonDown() {

}

void GlyphWidget::OnRightButtonUp() {

}

void GlyphWidget::Modified() {
    this->BuildRepresentation();
}
