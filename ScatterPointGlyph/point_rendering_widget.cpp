/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "point_rendering_widget.h"
#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include "vtkActor.h"
#include "vtkAssemblyNode.h"
#include "vtkAssemblyPath.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellPicker.h"
#include "vtkCommand.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPickingManager.h"
#include "vtkPlanes.h"
#include "vtkPointWidget.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkScalarBarActor.h"
#include "vtkLookupTable.h"
#include "vtkTextProperty.h"
#include "vtkPVScalarBarActor.h"
#include "scatter_point_dataset.h"

vtkStandardNewMacro(PointRenderingWidget);

PointRenderingWidget::PointRenderingWidget() {
	this->poly_data_ = vtkPolyData::New();
	this->mapper_ = vtkPolyDataMapper::New();
	this->actor_ = vtkActor::New();
	this->dataset_ = NULL;

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
	actor_->GetProperty()->SetPointSize(6);
	actor_->GetProperty()->SetColor(0.8, 0.8, 0.8);

    bar_actor_ = vtkPVScalarBarActor::New();
    //bar_actor_->SetNumberOfLabels(2);
    bar_actor_->SetMaximumWidthInPixels(20);
    bar_actor_->SetWidth(0.05);
    bar_actor_->GetLabelTextProperty()->SetFontFamilyToArial();
    bar_actor_->GetLabelTextProperty()->SetShadow(false);
    bar_actor_->GetLabelTextProperty()->SetBold(false);
    bar_actor_->GetLabelTextProperty()->SetColor(0.0, 0.0, 0.0);
    bar_actor_->GetLabelTextProperty()->SetFontSize(5);
    bar_actor_->GetTitleTextProperty()->SetFontFamilyToArial();
    bar_actor_->GetTitleTextProperty()->SetBold(true);
    bar_actor_->GetTitleTextProperty()->SetShadow(false);
    bar_actor_->GetTitleTextProperty()->SetColor(0.0, 0.0, 0.0);
    bar_actor_->GetTitleTextProperty()->SetFontSize(5);
    bar_actor_->SetPosition(0.91, 0.1);
    bar_actor_->SetHeight(0.7);
    bar_actor_->SetVisibility(false);

    scalar_lookup_table_ = vtkLookupTable::New();
    scalar_lookup_table_->SetValueRange(0, 1.0);
    scalar_lookup_table_->SetHueRange(0.8, 0.8);
    scalar_lookup_table_->SetAlphaRange(0.2, 1.0);
    scalar_lookup_table_->Build();
}

PointRenderingWidget::~PointRenderingWidget() {

}

void PointRenderingWidget::SetData(ScatterPointDataset* data) {
	dataset_ = data;

    this->BuildPointRepresentation();
}

void PointRenderingWidget::SetViewpot(float left, float right, float bottom, float top) {
    int point_count = 0;
    for (int i = 0; i < dataset_->point_num; ++i)
        if ((dataset_->original_point_pos[0][i] - left) * (dataset_->original_point_pos[0][i] - right) < 0
            && (dataset_->original_point_pos[1][i] - top) * (dataset_->original_point_pos[1][i] - bottom) < 0)
            point_count++;

    view_left_ = left;
    view_right_ = right;
    view_top_ = top;
    view_bottom_ = bottom;

    if (false)
        this->BuildDensityRepresentation();
    else
        this->BuildPointRepresentation();
}

void PointRenderingWidget::BuildPointRepresentation() {
    vtkSmartPointer< vtkPoints > original_points = vtkSmartPointer< vtkPoints >::New();
	vtkSmartPointer< vtkPolyData > original_point_poly = vtkSmartPointer< vtkPolyData >::New();
	original_point_poly->SetPoints(original_points);
	for (int i = 0; i < dataset_->original_point_pos[0].size(); ++i){
        if (!dataset_->is_valid[i]) continue;
		double new_pos[3];
		new_pos[0] = dataset_->original_point_pos[0][i];
		new_pos[1] = dataset_->original_point_pos[1][i];
		new_pos[2] = -0.1;
		original_points->InsertNextPoint(new_pos);
	}

	vtkSmartPointer< vtkVertexGlyphFilter > original_vertex_filter = vtkSmartPointer< vtkVertexGlyphFilter >::New();
	original_vertex_filter->SetInputData(original_point_poly);
	original_vertex_filter->Update();
	vtkSmartPointer< vtkPolyData > original_vertex_data = vtkSmartPointer< vtkPolyData >::New();

	poly_data_->Initialize();
	poly_data_->DeepCopy(original_vertex_filter->GetOutput());

	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
    colors->SetName("color");
	colors->SetNumberOfComponents(4);
	for (int i = 0; i < dataset_->original_point_pos[0].size(); ++i) {
        if (!dataset_->is_valid[i]) continue;
		colors->InsertNextTuple4(200, 200, 200, 150);
	}
	poly_data_->GetPointData()->SetScalars(colors);

	actor_->Modified();
}

void PointRenderingWidget::BuildDensityRepresentation() {
    
}

void PointRenderingWidget::SetEnabled(int enabling) {
	if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

	if (enabling) {
		if (this->Enabled) return;

        if (!this->CurrentRenderer) {
            this->SetCurrentRenderer(this->Interactor->FindPokedRenderer(
            this->Interactor->GetLastEventPosition()[0],
            this->Interactor->GetLastEventPosition()[1]));
            if (this->CurrentRenderer == NULL) return;
        }

		this->Enabled = 1;

		this->CurrentRenderer->AddActor(this->actor_);
        this->CurrentRenderer->AddActor(this->bar_actor_);

        this->InvokeEvent(vtkCommand::EnableEvent,NULL);
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->CurrentRenderer->RemoveActor(this->actor_);
        this->CurrentRenderer->RemoveActor(this->bar_actor_);

        this->InvokeEvent(vtkCommand::DisableEvent,NULL);
        this->SetCurrentRenderer(NULL);
	}

	this->Interactor->Render();
}

void PointRenderingWidget::SetSelectedIds(vector<vector<int>>& ids) {
    vtkUnsignedCharArray* colors = (vtkUnsignedCharArray*)poly_data_->GetPointData()->GetScalars("color");
    // clear previous selection
    for (int i = 0; i < selected_ids_.size(); ++i) 
        for (int j = 0; j < selected_ids_[i].size(); ++j) {
            colors->SetTuple4(selected_ids_[i][j], 200, 200, 200, 255);
	    }

    // add new selection
    for (int i = 0; i < ids.size(); ++i) 
        for (int j = 0; j < ids[i].size(); ++j) {
            colors->SetTuple4(ids[i][j], 255, 140, 61, 255);
	    }

    selected_ids_ = ids;

    colors->Modified();
    this->Interactor->Render();
}

void PointRenderingWidget::SetColorMappingOn(QString name, vector<float>& vals, float min_val, float max_val) {
    scalar_lookup_table_->SetValueRange(min_val, max_val);
    scalar_lookup_table_->SetTableRange(min_val, max_val);
    scalar_lookup_table_->SetHueRange(0.0, 0.0);
    scalar_lookup_table_->SetSaturationRange(0.0, 0.0);
    scalar_lookup_table_->SetValueRange(0.7, 0.3);
    scalar_lookup_table_->Build();

    bar_actor_->SetLookupTable(scalar_lookup_table_);
    bar_actor_->SetVisibility(true);
    bar_actor_->SetTitle(name.toLocal8Bit().data());
    bar_actor_->Modified();

    vtkUnsignedCharArray* color_array = (vtkUnsignedCharArray*)poly_data_->GetPointData()->GetScalars("color");
    double rgb[3];
    int test_num = color_array->GetNumberOfTuples();
	for (int i = 0; i < vals.size() && i < color_array->GetNumberOfTuples(); ++i) {
        scalar_lookup_table_->GetColor(vals[i], rgb);
		color_array->SetTuple4(i, rgb[0] * 255, rgb[1] * 255, rgb[2] * 255, 255);
	}
    color_array->Modified();
    this->Interactor->Render();
}

void PointRenderingWidget::SetColorMappingOff() {
    vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
    double rgb[3];
	for (int i = 0; i < color_array->GetNumberOfTuples(); ++i) {
		color_array->SetTuple4(i, 200, 200, 200, 255);
	}
    color_array->Modified();

    bar_actor_->SetVisibility(false);
    bar_actor_->Modified();

    this->Interactor->Render();
}