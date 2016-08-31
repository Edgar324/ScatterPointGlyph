/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "map_rendering_widget.h"
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

vtkStandardNewMacro(MapRenderingWidget);

MapRenderingWidget::MapRenderingWidget() {
    this->actor_ = vtkActor::New();
	this->polydata_ = vtkPolyData::New();
	this->data_mapper_ = vtkPolyDataMapper::New();
	this->data_mapper_->SetInputData(this->polydata_);
	this->actor_->SetMapper(this->data_mapper_);
    this->actor_->GetProperty()->SetColor(0.5, 0.5, 0.5);
}

MapRenderingWidget::~MapRenderingWidget() {

}

void MapRenderingWidget::SetEnabled(int enabled) {
    if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

	if (enabled) {
		if (this->Enabled) return;

        if (!is_data_loaded_) this->BuildRepresentation();

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

void MapRenderingWidget::BuildRepresentation() {
    polydata_->Initialize();
    vtkPoints* points = vtkPoints::New();
    polydata_->SetPoints(points);
    vtkCellArray* poly_array = vtkCellArray::New();
	polydata_->SetLines(poly_array);

    std::ifstream border_file("./Resources/border.txt");

    if ( border_file.good() ){
        while ( !border_file.eof() ){
            int poly_size;
            float x, y;
            std::vector< float > temp_poly;
            border_file >> poly_size;
            temp_poly.resize(poly_size * 2);
            bool is_negative = false;
            for ( int i = 0; i < poly_size; ++i ){
                border_file >> temp_poly[2 * i] >> temp_poly[2 * i + 1];
                if ( temp_poly[2 * i] < 0 ) is_negative = true;
            }
            /*if ( is_negative )
                for ( int i = 0; i < temp_poly.size() / 2; ++i ) temp_poly[2 * i] += 360;*/
            bool is_used = false;
            /*for ( int i = 0; i < temp_poly.size() / 2; ++i )
                if ( (temp_poly[2 * i] - start_x) * (temp_poly[2 * i] - end_x) <= 0 && (temp_poly[2 * i +1] - start_y) * (temp_poly[2 * i + 1] - end_y) <= 0 ){
                    is_used = true;
                    break;
                }*/
            is_used = true;
            std::vector< vtkIdType > poly_ids;
            if (is_used) {
                for (int i = 0; i < temp_poly.size() / 2; ++i) {
                    poly_ids.push_back(points->InsertNextPoint(temp_poly[i * 2], temp_poly[i * 2 + 1], 0));
                }
                polydata_->InsertNextCell(VTK_POLY_LINE, poly_ids.size(), poly_ids.data());
            }
        }
    }

    polydata_->Modified();
    actor_->Modified();

    is_data_loaded_ = true;
}