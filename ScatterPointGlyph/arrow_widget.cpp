/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "arrow_widget.h"
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

vtkStandardNewMacro(ArrowWidget);

ArrowWidget::ArrowWidget() {
    this->arrow_actor_ = vtkActor::New();
	this->arrow_poly_ = vtkPolyData::New();
	this->arrow_data_mapper_ = vtkPolyDataMapper::New();
	this->arrow_data_mapper_->SetInputData(this->arrow_poly_);
	this->arrow_actor_->SetMapper(this->arrow_data_mapper_);
    this->arrow_actor_->GetProperty()->SetColor(1.0, 0.55, 0.24);
}

ArrowWidget::~ArrowWidget() {

}

void ArrowWidget::SetData(vector<float>& points, float node_radius) {
    point_pos_ = points;
    node_radius_ = node_radius;

    this->BuildRepresentation();
}

void ArrowWidget::SetEnabled(int enabled) {
    if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

	if (enabled) {
		if (this->Enabled) return;

        if (!this->CurrentRenderer) {
            this->SetCurrentRenderer(this->Interactor->FindPokedRenderer(
            this->Interactor->GetLastEventPosition()[0],
            this->Interactor->GetLastEventPosition()[1]));
            if (this->CurrentRenderer == NULL) return;
        }

		this->Enabled = 1;

		this->CurrentRenderer->AddActor(this->arrow_actor_);

        this->InvokeEvent(vtkCommand::EnableEvent,NULL);
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->CurrentRenderer->RemoveActor(this->arrow_actor_);

        this->InvokeEvent(vtkCommand::DisableEvent,NULL);
        this->SetCurrentRenderer(NULL);
	}

	this->Interactor->Render();
}

void ArrowWidget::BuildRepresentation() {
    this->arrow_poly_->Initialize();

	vtkPoints* linked_points = vtkPoints::New();
	arrow_poly_->SetPoints(linked_points);

	vtkCellArray* linked_poly_array = vtkCellArray::New();
	arrow_poly_->SetPolys(linked_poly_array);

    vtkIdType cell_ids[4];

	for (int i = 0; i < point_pos_.size() / 4; ++i){
		float dir_vec[2], orth_vec[2];
		dir_vec[0] = point_pos_[4 * i + 2] - point_pos_[4 * i];
		dir_vec[1] = point_pos_[4 * i + 3] - point_pos_[4 * i + 1];
		float dir_length = sqrt(pow(dir_vec[0], 2) + pow(dir_vec[1], 2));
		dir_vec[0] /= dir_length;
		dir_vec[1] /= dir_length;
		orth_vec[0] = -1 * dir_vec[1];
		orth_vec[1] = dir_vec[0];

		float begin_pos[2], end_pos[2], arrow_end_pos[2];
		begin_pos[0] = point_pos_[4 * i] + node_radius_ * dir_vec[0] * 1.2;
		begin_pos[1] = point_pos_[4 * i + 1] + node_radius_ * dir_vec[1] * 1.2;
		end_pos[0] = point_pos_[4 * i + 2] - node_radius_ * dir_vec[0] * 1.5;
        end_pos[1] = point_pos_[4 * i + 3] - node_radius_ * dir_vec[1] * 1.5;
        arrow_end_pos[0] = point_pos_[4 * i + 2] - node_radius_ * dir_vec[0] * 1.2;
        arrow_end_pos[1] = point_pos_[4 * i + 3] - node_radius_ * dir_vec[1] * 1.2;

		float value_dis = 0.05;

		cell_ids[0] = linked_points->InsertNextPoint(begin_pos[0] + node_radius_ * orth_vec[0] * value_dis, begin_pos[1] + node_radius_ * orth_vec[1] * value_dis, 0.001);
		cell_ids[1] = linked_points->InsertNextPoint(end_pos[0] + node_radius_ * orth_vec[0] * value_dis, end_pos[1] + node_radius_ * orth_vec[1] * value_dis, 0.001);
		cell_ids[2] = linked_points->InsertNextPoint(end_pos[0] - node_radius_ * orth_vec[0] * value_dis, end_pos[1] - node_radius_ * orth_vec[1] * value_dis, 0.001);
		cell_ids[3] = linked_points->InsertNextPoint(begin_pos[0] - node_radius_ * orth_vec[0] * value_dis, begin_pos[1] - node_radius_ * orth_vec[1] * value_dis, 0.001);
		arrow_poly_->InsertNextCell(VTK_POLYGON, 4, cell_ids);

        value_dis = 0.2;
        cell_ids[0] = linked_points->InsertNextPoint(end_pos[0] + node_radius_ * orth_vec[0] * value_dis, end_pos[1] + node_radius_ * orth_vec[1] * value_dis, 0.001);
        cell_ids[1] = linked_points->InsertNextPoint(arrow_end_pos[0], arrow_end_pos[1], 0.001);
		cell_ids[2] = linked_points->InsertNextPoint(end_pos[0] - node_radius_ * orth_vec[0] * value_dis, end_pos[1] - node_radius_ * orth_vec[1] * value_dis, 0.001);
        arrow_poly_->InsertNextCell(VTK_POLYGON, 3, cell_ids);
	}

	this->arrow_actor_->Modified();
}

