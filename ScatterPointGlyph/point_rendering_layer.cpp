#include "point_rendering_layer.h"
#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>

PointRenderingLayer::PointRenderingLayer() 
	: BasicGlyphLayer() {

}

PointRenderingLayer::~PointRenderingLayer() {

}

void PointRenderingLayer::SetData(vtkPolyData* data) {
	poly_data_ = data;

	UpdatePointGlyph();
}

void PointRenderingLayer::SetHighlightPointIndex(std::vector< int >& index) {
	point_index_ = index;

	vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
	for (int i = 0; i < poly_data_->GetPoints()->GetNumberOfPoints(); ++i)
		color_array->SetTuple3(i, 128, 128, 128);
	for (int i = 0; i < index.size(); ++i) color_array->SetTuple3(index[i], 128, 0, 0);
}

void PointRenderingLayer::UpdatePointGlyph() {
	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	for (int i = 0; i < poly_data_->GetPoints()->GetNumberOfPoints(); ++i)
		colors->InsertNextTuple3(128, 128, 128);

	poly_data_->GetPointData()->SetScalars(colors);

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
	actor_->GetProperty()->SetPointSize(3);
}