#include "field_rendering_layer.h"
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
#include <vtkLookupTable.h>
#include <vtkDelaunay2D.h>
#include "scatter_point_dataset.h"

FieldRenderingLayer::FieldRenderingLayer() {
	this->is_segment_probability_ = false;
	this->dataset_ = NULL;
	this->var_index_ = -1;

	this->point_polydata_ = vtkPolyData::New();
	this->point_mapper_ = vtkPolyDataMapper::New();
	this->point_actor_ = vtkActor::New();
	this->point_mapper_->SetInputData(this->point_polydata_);
	this->point_actor_->SetMapper(this->point_mapper_);
	this->point_actor_->GetProperty()->SetPointSize(6);

	this->field_polydata_ = vtkPolyData::New();
	this->field_poly_mapper_ = vtkPolyDataMapper::New();
	this->field_actor_ = vtkActor::New();
	this->field_poly_mapper_->SetInputData(this->field_polydata_);
	this->field_actor_->SetMapper(this->field_poly_mapper_);

	color_bar_actor_ = vtkScalarBarActor::New();
	color_bar_actor_->SetMaximumWidthInPixels(70);
	lookup_table_ = vtkLookupTable::New();
	lookup_table_->SetTableRange(0, 10);
	this->color_bar_actor_->SetNumberOfLabels(2);
	this->color_bar_actor_->SetLookupTable(lookup_table_);
}

FieldRenderingLayer::~FieldRenderingLayer() {

}

void FieldRenderingLayer::SetFieldData(ScatterPointDataset* data, int var_index) {
	this->dataset_ = data;
	this->var_index_ = var_index;
	this->is_segment_probability_ = false;

	ConstructFieldActors();
}

void FieldRenderingLayer::SetColorScalar(std::vector< float >& values, std::vector< double >& rgb) {
	lookup_table_->SetTableRange(values[0], values[values.size() - 1]);
	lookup_table_->SetNumberOfTableValues(values.size() - 1);
	for (int i = 0; i < values.size() - 1; ++i) {
		lookup_table_->SetTableValue(i, &rgb[4 * i]);
	}
	lookup_table_->Modified();
	this->color_bar_actor_->SetLookupTable(lookup_table_);
	color_bar_actor_->SetNumberOfLabels(values.size());
	color_bar_actor_->Modified();

	this->point_mapper_->ScalarVisibilityOn();
	this->point_mapper_->SetScalarModeToUsePointData();
	this->point_mapper_->SetColorModeToMapScalars();
	this->point_mapper_->UseLookupTableScalarRangeOn();
	this->point_mapper_->SetLookupTable(lookup_table_);
	this->point_actor_->Modified();

	this->field_poly_mapper_->ScalarVisibilityOn();
	this->field_poly_mapper_->SetScalarModeToUsePointData();
	this->field_poly_mapper_->SetColorModeToMapScalars();
	this->field_poly_mapper_->UseLookupTableScalarRangeOn();
	this->field_poly_mapper_->SetLookupTable(lookup_table_);
	this->field_actor_->Modified();
}

void FieldRenderingLayer::ConstructFieldActors() {
	// construct point field
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

	this->point_polydata_->Initialize();
	this->point_polydata_->DeepCopy(original_vertex_filter->GetOutput());

	vtkFloatArray* value_array = vtkFloatArray::New();
	value_array->SetNumberOfComponents(1);
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i) {
		value_array->InsertNextTuple1(dataset_->original_point_values[i][var_index_]);
	}
	this->point_polydata_->GetPointData()->SetScalars(value_array);
	this->point_actor_->Modified();

	// construct surface field
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer< vtkDelaunay2D >::New();
	delaunay->SetInputData(original_point_poly);
	delaunay->SetTolerance(0.00001);
	delaunay->SetBoundingTriangulation(false);
	delaunay->Update();
	this->field_polydata_->DeepCopy(delaunay->GetOutput());
	this->field_polydata_->GetPointData()->SetScalars(value_array);
	this->field_actor_->Modified();
}

void FieldRenderingLayer::SetEnabled(int enabling) {
	if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

	if (enabling) {
		if (this->Enabled) return;

		if (!this->CurrentRenderer) {
			this->SetCurrentRenderer(
				this->Interactor->FindPokedRenderer(
				this->Interactor->GetLastEventPosition()[0],
				this->Interactor->GetLastEventPosition()[1]));
			if (this->CurrentRenderer == NULL) return;
		}

		this->Enabled = 1;

		this->CurrentRenderer->AddActor(this->point_actor_);
		this->CurrentRenderer->AddActor(this->color_bar_actor_);
		this->CurrentRenderer->AddActor(this->field_actor_);
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->CurrentRenderer->RemoveActor(this->point_actor_);
		this->CurrentRenderer->RemoveActor(this->color_bar_actor_);
		this->CurrentRenderer->RemoveActor(this->field_actor_);

		this->InvokeEvent(vtkCommand::DisableEvent, NULL);
		this->SetCurrentRenderer(NULL);
	}

	this->Interactor->Render();
}