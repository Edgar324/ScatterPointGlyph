/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "density_rendering_widget.h"
#include <algorithm>
using namespace std;

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
#include "color_mapping_generator.h"

vtkStandardNewMacro(DensityRenderingWidget);

DensityRenderingWidget::DensityRenderingWidget() {
    this->actor_ = vtkActor::New();
	this->polydata_ = vtkPolyData::New();
	this->data_mapper_ = vtkPolyDataMapper::New();
	this->data_mapper_->SetInputData(this->polydata_);
	this->actor_->SetMapper(this->data_mapper_);
    this->actor_->GetProperty()->SetColor(1.0, 0.55, 0.24);
}

DensityRenderingWidget::~DensityRenderingWidget() {

}

void DensityRenderingWidget::SetData(vector<vector<float>>& pos, vector<int>& cluster_index) {
    point_pos_ = pos;
    cluster_index_ = cluster_index;

    this->BuildRepresentation();

    this->Interactor->Render();
}

void DensityRenderingWidget::SetEnabled(int enabled) {
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

		this->CurrentRenderer->AddActor(this->actor_);

        this->InvokeEvent(vtkCommand::EnableEvent,NULL);
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->CurrentRenderer->RemoveActor(this->actor_);

        this->InvokeEvent(vtkCommand::DisableEvent,NULL);
        this->SetCurrentRenderer(NULL);
	}

	this->Interactor->Render();
}

void DensityRenderingWidget::BuildRepresentation() {
    // construct adaptive rate
    int point_num = point_pos_[0].size();

	vector<vector<float>> point_dis;
	point_dis.resize(point_num);
	for (int i = 0; i < point_num; ++i)
		point_dis[i].resize(point_num);
	for (int i = 0; i < point_num - 1; ++i) {
		point_dis[i][i] = 0;
		for (int j = i + 1; j < point_num; ++j) {
			float dis = sqrt(pow(point_pos_[0][i] - point_pos_[0][j], 2) + pow(point_pos_[1][i] - point_pos_[1][j], 2));
			point_dis[i][j] = dis;
			point_dis[j][i] = dis;
		}
	}
    float max_dis = 0;
    for (int i = 0; i < point_dis.size(); ++i) {
        sort(point_dis[i].begin(), point_dis[i].end());
        if (point_dis[i][point_dis[i].size() - 1] > max_dis)
            max_dis = point_dis[i][point_dis[i].size() - 1];
    }

	adaptive_rate_.resize(point_num);
	int kNum = 20;
	for (int i = 0; i < point_num; ++i) {
		float ave_dis = 0;
		for (int j = 0; j < kNum; ++j)
			ave_dis += point_dis[i][j];
		ave_dis /= (kNum - 1);
		if (ave_dis > 0.1 * max_dis) ave_dis = 0.1 * max_dis;
		adaptive_rate_[i] = ave_dis;
	}

    vector<QColor> mapping_colors;
    ColorMappingGenerator::GetInstance()->GetQualitativeColors(10, mapping_colors);

    this->polydata_->Initialize();

    vtkPoints* points = vtkPoints::New();
	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
	colors->SetNumberOfComponents(4);

	polydata_->Initialize();
	polydata_->SetPoints(points);
	polydata_->GetPointData()->SetScalars(colors);

	vtkCellArray* poly_array = vtkCellArray::New();
	polydata_->SetPolys(poly_array);

    vtkIdType ids[21];
    for (int i = 0; i < point_num; ++i) {
        float x = point_pos_[0][i];
        float y = point_pos_[1][i];
        int cluster = cluster_index_[i];
        if (cluster > 10) cluster = 10;

        ids[0] = points->InsertNextPoint(x, y, -0.0001);
        colors->InsertNextTuple4(mapping_colors[cluster].red(), mapping_colors[cluster].green(), mapping_colors[cluster].blue(), 80);
      
		float radius = adaptive_rate_[i];
        for (int j = 0; j < 20; ++j) {
            float end_arc = j * 2 * 3.14159 / 19;
            ids[j + 1] = points->InsertNextPoint(x + radius * cos(end_arc), y + radius * sin(end_arc), -0.0001);
            colors->InsertNextTuple4(mapping_colors[cluster].red(), mapping_colors[cluster].green(), mapping_colors[cluster].blue(), 10);
        }

        vtkIdType poly_ids[3];
        for (int j = 0; j < 19; ++j) {
            poly_ids[0] = ids[0];
            poly_ids[1] = ids[j + 1];
            poly_ids[2] = ids[j + 2];

            polydata_->InsertNextCell(VTK_TRIANGLE, 3, poly_ids);
        }
    }

    polydata_->Modified();
    actor_->Modified();
}