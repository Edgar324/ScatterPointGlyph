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

FieldRenderingLayer::FieldRenderingLayer() {
	this->is_segment_probability_ = false;
	this->dataset_ = NULL;
	this->var_index_ = -1;

	this->point_polydata_ = vtkPolyData::New();
	this->point_mapper_ = vtkPolyDataMapper::New();
	this->point_actor_ = vtkActor::New();
	this->point_mapper_->SetInputData(this->point_polydata_);
	this->point_actor_->SetMapper(this->point_mapper_);

	this->field_polydata_ = vtkPolyData::New();
	this->field_poly_mapper_ = vtkPolyDataMapper::New();
	this->field_actor_ = vtkActor::New();
	this->field_poly_mapper_->SetInputData(this->field_polydata_);
	this->field_actor_->SetMapper(this->field_poly_mapper_);
}

FieldRenderingLayer::~FieldRenderingLayer() {

}

void FieldRenderingLayer::SetFieldData(ScatterPointDataset* data, int var_index) {
	this->dataset_ = data;
	this->var_index_ = var_index;
	this->is_segment_probability_ = false;
}

void FieldRenderingLayer::ConstructFieldActors() {
	
}

void FieldRenderingLayer::SetEnabled(int enabling) {

}