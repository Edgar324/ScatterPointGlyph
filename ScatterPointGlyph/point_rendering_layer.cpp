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

void PointRenderingLayer::SetPointValue(std::vector< float >& values){
	vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
	for (int i = 0; i < values.size(); ++i) {
		int grey = (int)(values[i] * 255);
		color_array->SetTuple3(i, grey, grey, grey);
	}
	color_array->Modified();
}

void PointRenderingLayer::SetHighlightPointIndex(std::vector< int >& index) {
	point_index_ = index;

	if (index.size() <= 2) return;

	int red = 255 * rand() / RAND_MAX;
	int green = 255 * rand() / RAND_MAX;
	int blue = 255 * rand() / RAND_MAX;
	vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
	for (int i = 0; i < index.size(); ++i) color_array->SetTuple3(index[i], red, green, blue);
	color_array->Modified();
}

void PointRenderingLayer::SetClusterIndex(int cluster_count, std::vector< int >& point_index) {
	std::vector< int > rgb;
	rgb.resize(3 * cluster_count);
	for (int i = 0; i < cluster_count; ++i) {
		
		rgb[3 * i] = 255 * rand() / RAND_MAX;
		rgb[3 * i + 1] = 255 * rand() / RAND_MAX;
		rgb[3 * i + 2] = 255 * rand() / RAND_MAX;
	}

	vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());
	for (int i = 0; i < point_index.size(); ++i) {
		int cindex = point_index[i] * 3;
		if (cindex < 0)
			color_array->SetTuple3(i, 0, 0, 0);
		else
			color_array->SetTuple3(i, rgb[cindex], rgb[cindex + 1], rgb[cindex + 2]);
	}
	color_array->Modified();
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
	actor_->GetProperty()->SetPointSize(6);

	actor_->Modified();
}