#include "map_rendering_layer.h"
#include <vtkContourFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPlaneSource.h>
#include <vtkMarchingContourFilter.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkStripper.h>
#include <vtkCellArray.h>
#include <vtkContourTriangulator.h>
#include <vtkAppendPolyData.h>
#include <vtkProperty.h>
#include <vtkCleanPolyData.h>

MapRenderingLayer::MapRenderingLayer() 
	: BasicGlyphLayer() {
	
}

MapRenderingLayer::~MapRenderingLayer() {

}

void MapRenderingLayer::SetFieldData(float* data, int w, int h, float wSpacing, float hSpacing, float startX, float startY) {
	vtkSmartPointer< vtkMarchingContourFilter > contours = vtkSmartPointer< vtkMarchingContourFilter >::New();
	vtkSmartPointer< vtkPlaneSource > plane = vtkSmartPointer< vtkPlaneSource >::New();
	plane->SetXResolution(w + 1);
	plane->SetYResolution(h + 1);
	plane->SetCenter(0, 0, 0);
	plane->SetNormal(0, 0, 1);
	plane->Update();

	vtkSmartPointer< vtkFloatArray > map_scalars = vtkSmartPointer< vtkFloatArray >::New();
	map_scalars->SetNumberOfComponents(1);
	int temp_index = 0;
	for (int i = 0; i <= h + 1; ++i)
		for (int j = 0; j <= w + 1; ++j){
			if (i == 0 || i == h + 1
				|| j == 0 || j == w + 1)
				map_scalars->InsertNextTuple1(-1e10);
			else{
				map_scalars->InsertNextTuple1(*(data + temp_index));
				temp_index++;
			}
		}
	plane->GetOutput()->GetPointData()->SetScalars(map_scalars);

	if (iso_contours_.size() != 0) {
		for (int i = 0; i < iso_contours_.size(); ++i) iso_contours_[i]->Delete();
		iso_contours_.clear();
	}

	std::vector< float > iso_values;
	for (int i = 0; i < 5; ++i) iso_values.push_back(2 * (i + 1));

	for (int i = 0; i < iso_values.size(); ++i){
		contours->SetInputConnection(plane->GetOutputPort());
		contours->GenerateValues(1, iso_values[i], 9999);
		contours->Update();

		vtkPolyData* contour_poly = contours->GetOutput();
		vtkPolyData* new_data = vtkPolyData::New();
		new_data->DeepCopy(contour_poly);

		for (int j = 0; j < contour_poly->GetNumberOfPoints(); ++j) {
			double* p = contour_poly->GetPoint(j);
			double new_pos[3];
			float x = p[0] * (float)(w + 1) / (w - 1);
			float y = p[1] * (float)(h + 1) / (h - 1);
			new_pos[0] = (x + 0.5) * (w - 1) * wSpacing + startX;
			new_pos[1] = (y + 0.5) * (h - 1) * hSpacing + startY;
			new_pos[2] = 0.1;
			new_data->GetPoints()->SetPoint(j, new_pos);
		}

		iso_contours_.push_back(new_data);
	}

	this->Update();
}

void MapRenderingLayer::Update() {
	if (poly_data_ != NULL) {
		poly_data_->Delete();
	}
	poly_data_ = vtkPolyData::New();

	vtkSmartPointer< vtkAppendPolyData > appendFilter =
		vtkSmartPointer< vtkAppendPolyData >::New();
	//appendFilter->AddInputData(border_data_);
	for (int i = 0; i < iso_contours_.size(); ++i) appendFilter->AddInputData(iso_contours_[i]);
	vtkSmartPointer< vtkCleanPolyData > cleanFilter =
		vtkSmartPointer< vtkCleanPolyData >::New();
	cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
	cleanFilter->Update();
	poly_data_->DeepCopy(cleanFilter->GetOutput());

	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");
	for (int i = 0; i < poly_data_->GetPoints()->GetNumberOfPoints(); ++i)
		colors->InsertNextTuple3(255, 0, 0);
	poly_data_->GetPointData()->SetScalars(colors);

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
	actor_->GetProperty()->SetLineWidth(2.0);
}