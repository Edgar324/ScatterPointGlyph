/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "indicator_widget.h"

#include <vtkObjectFactory.h>
#include <vtk3DWidget.h>
#include <vtkCommand.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellPicker.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
#include <vtkTextProperty.h>

vtkStandardNewMacro(IndicatorWidget);

IndicatorWidget::IndicatorWidget() {
    this->actor_ = vtkActor::New();
	this->poly_ = vtkPolyData::New();
	this->data_mapper_ = vtkPolyDataMapper::New();
	this->data_mapper_->SetInputData(this->poly_);
	this->actor_->SetMapper(this->data_mapper_);
    this->actor_->GetProperty()->SetColor(0.5, 0.5, 0.5);

    indicator_text_ = vtkTextActor3D::New();
    indicator_text_->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
	indicator_text_->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
	indicator_text_->GetTextProperty()->SetBold(true);
	indicator_text_->GetTextProperty()->SetFontFamilyToArial();
}

IndicatorWidget::~IndicatorWidget() {

}

void IndicatorWidget::SetData(int max_count) {
    this->max_count_ = max_count;

    this->BuildRepresentation();
}

void IndicatorWidget::SetRenderer(vtkRenderer* renderer) {
    this->renderer_ = renderer;
}


void IndicatorWidget::SetEnabled(int enabled) {
    if (!this->renderer_) {
		vtkErrorMacro(<< "The renderer must be set prior to enabling/disabling widget");
		return;
	}

	if (enabled) {
		if (this->Enabled) return;

        this->SetCurrentRenderer(this->renderer_);

		this->Enabled = 1;

		this->CurrentRenderer->AddActor(this->actor_);
        this->CurrentRenderer->AddActor(indicator_text_);

        this->CurrentRenderer->ResetCamera();

        this->InvokeEvent(vtkCommand::EnableEvent,NULL);
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->CurrentRenderer->RemoveActor(this->actor_);
        this->CurrentRenderer->RemoveActor(indicator_text_);

        this->InvokeEvent(vtkCommand::DisableEvent,NULL);
        this->SetCurrentRenderer(NULL);
	}

	this->Interactor->Render();
}

void IndicatorWidget::BuildRepresentation() {
    vtkPoints* indicator_points = vtkPoints::New();
	poly_->Initialize();
	poly_->SetPoints(indicator_points);
	vtkCellArray* indicator_poly_array = vtkCellArray::New();
	poly_->SetPolys(indicator_poly_array);

    int seg_per_circle = 50;
    int seg_num = seg_per_circle - 1;
    float gray_value = 128;
    std::vector<vtkIdType > inner_ids, outer_ids;
    for (int j = 0; j <= seg_num; ++j) {
        float end_arc = -1 * j * 3.14159 * 2 / seg_per_circle + 3.14159 * 0.5;
        float x = 1.0 * cos(end_arc) * 1.05;
		float y = 1.0 * sin(end_arc) * 1.05;

        vtkIdType id_one = indicator_points->InsertNextPoint(x, y, 0.000);
        inner_ids.push_back(id_one);

        x = 1.0 * cos(end_arc) * 1.18;
		y = 1.0 * sin(end_arc) * 1.18;

        vtkIdType id_two = indicator_points->InsertNextPoint(x, y, 0.000);
        outer_ids.push_back(id_two);
    }

    vtkIdType cell_ids[5];
    for (int j = 0; j < outer_ids.size() - 1; ++j) {
        cell_ids[0] = outer_ids[j];
        cell_ids[1] = outer_ids[j + 1];
        cell_ids[2] = inner_ids[j];
        poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

        cell_ids[0] = outer_ids[j + 1];
        cell_ids[1] = inner_ids[j + 1];
        cell_ids[2] = inner_ids[j];
        poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
    }
    actor_->Modified();

    char buffer[20];
	itoa(this->max_count_, buffer, 10);
    std::string str = std::string("Max Point Num: ") + std::string(buffer);
	indicator_text_->SetInput(str.c_str());
	indicator_text_->SetPosition(-1, -1.5, 0);
	indicator_text_->GetTextProperty()->SetFontSize(20);
	float scale = 0.2 / 20;
	indicator_text_->SetScale(scale, scale, scale);
	indicator_text_->Modified();
}
