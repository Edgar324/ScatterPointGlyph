/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "glyph_widget.h"

vtkStandardNewMacro(GlyphWidget);

GlyphWidget::GlyphWidget() {

}

GlyphWidget::~GlyphWidget() {

}

void GlyphWidget::SetData(int cluster_index, vector<QString>& names, 
    vector<float>& means, vector<float>& std_devs, 
    float saliency, int node_count, int max_node_count) {
    this->cluster_index_ = cluster_index;
    this->names_ = names;
    this->means_ = means;
    this->std_devs_ = std_devs;
    this->saliency_ = saliency;
    this->node_count_ = node_count;
    this->max_node_count_ = max_node_count;
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
    
}

void GlyphWidget::OnLeftButtonDown() {

}

void GlyphWidget::OnLeftButtonUp() {

}

void GlyphWidget::OnRightButtonDown() {

}

void GlyphWidget::OnRightButtonUp() {

}
