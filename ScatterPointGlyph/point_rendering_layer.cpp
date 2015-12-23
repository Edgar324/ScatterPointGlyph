#include "point_rendering_layer.h"
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
#include "scatter_point_dataset.h"

PointRenderingLayer::PointRenderingLayer() {
	this->poly_data_ = vtkPolyData::New();
	this->mapper_ = vtkPolyDataMapper::New();
	this->actor_ = vtkActor::New();

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
	actor_->GetProperty()->SetPointSize(6);
	actor_->GetProperty()->SetColor(0.5, 0.5, 0.5);

	is_category_on_ = false;
	this->current_selection_index_ = -1;
}

PointRenderingLayer::~PointRenderingLayer() {

}

void PointRenderingLayer::SetData(ScatterPointDataset* data) {
	dataset_ = data;

	vtkSmartPointer< vtkPoints > original_points = vtkSmartPointer< vtkPoints >::New();
	vtkSmartPointer< vtkPolyData > original_point_poly = vtkSmartPointer< vtkPolyData >::New();
	original_point_poly->SetPoints(original_points);
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i){
		double new_pos[3];
		new_pos[0] = dataset_->original_point_pos[i][0];
		new_pos[1] = dataset_->original_point_pos[i][1];
		new_pos[2] = 0;
		original_points->InsertNextPoint(new_pos);
	}

	vtkSmartPointer< vtkVertexGlyphFilter > original_vertex_filter = vtkSmartPointer< vtkVertexGlyphFilter >::New();
	original_vertex_filter->SetInputData(original_point_poly);
	original_vertex_filter->Update();
	vtkSmartPointer< vtkPolyData > original_vertex_data = vtkSmartPointer< vtkPolyData >::New();

	poly_data_->Initialize();
	poly_data_->DeepCopy(original_vertex_filter->GetOutput());

	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
	colors->SetNumberOfComponents(3);
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i) {
		colors->InsertNextTuple3(128, 128, 128);
	}
	poly_data_->GetPointData()->SetScalars(colors);

	actor_->Modified();
}

void PointRenderingLayer::SetEnabled(int enabling) {
	if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

	if (enabling) {
		if (this->Enabled) return;

		this->Enabled = 1;

		this->DefaultRenderer->AddActor(this->actor_);
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->DefaultRenderer->RemoveActor(this->actor_);

		this->InvokeEvent(vtkCommand::DisableEvent, NULL);
	}

	this->Interactor->Render();
}

void PointRenderingLayer::SetPointValue(std::vector< float >& values){
	point_values_ = values;

	vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
	for (int i = 0; i < values.size(); ++i) {
		int grey = (int)((1.0 - values[i]) * 255);
		color_array->SetTuple3(i, grey, grey, grey);
	}
	color_array->Modified();
}

void PointRenderingLayer::SetClusterIndex(int cluster_count, std::vector< int >& point_index) {
	this->cluster_count_ = cluster_count;
	this->point_index_ = point_index;

	srand((unsigned int)time(0));
	cluster_color_.resize(3 * cluster_count);
	for (int i = 0; i < cluster_count; ++i) {
		cluster_color_[3 * i] = 255 * rand() / RAND_MAX;
		cluster_color_[3 * i + 1] = 255 * rand() / RAND_MAX;
		cluster_color_[3 * i + 2] = 255 * rand() / RAND_MAX;
	}

	if (is_category_on_) SetCategoryOn();
}

void PointRenderingLayer::SetHighlightCluster(int index) {
	if (is_category_on_) SetCategoryOn();
	else SetCategoryOff();

	this->current_selection_index_ = index;

	if (this->current_selection_index_ != -1) {
		vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
		for (int i = 0; i < point_index_.size(); ++i)
			if (point_index_[i] == this->current_selection_index_) {
				color_array->SetTuple3(i, 255, 0.0, 0.0);
			}
		color_array->Modified();
	}
}

void PointRenderingLayer::SetCategoryOn() {
	is_category_on_ = true;

	vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
	for (int i = 0; i < point_index_.size(); ++i) {
		int cindex = point_index_[i] * 3;
		if (cindex < 0)
			color_array->SetTuple3(i, 0, 0, 0);
		else
			color_array->SetTuple3(i, cluster_color_[cindex], cluster_color_[cindex + 1], cluster_color_[cindex + 2]);
	}
	color_array->Modified();
}

void PointRenderingLayer::SetCategoryOff() {
	is_category_on_ = false;

	vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
	if (point_values_.size() != 0) {
		for (int i = 0; i < point_values_.size(); ++i) {
			int grey = (int)((1.0 - point_values_[i]) * 255);
			color_array->SetTuple3(i, grey, grey, grey);
		}
	} else {
		for (int i = 0; i < dataset_->original_point_pos.size(); ++i) {
			color_array->SetTuple3(i, 128, 128, 128);
		}
	}
	color_array->Modified();
}

std::vector< bool > PointRenderingLayer::GetPointSelection() {
	std::vector< bool > selection;
	if (this->current_selection_index_ != -1 && this->point_index_.size() != 0) {
		selection.resize(dataset_->original_point_pos.size(), false);
		for (int i = 0; i < point_index_.size(); ++i)
			if (point_index_[i] == this->current_selection_index_) selection[i] = true;
	}
	return selection;
}