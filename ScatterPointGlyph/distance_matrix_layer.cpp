#include "distance_matrix_layer.h"
#include "vtkTextActor.h"
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
#include "vtkProperty2D.h"
#include "vtkTextWidget.h"
#include "vtkTextRepresentation.h"
#include "vtkTextActor3D.h"
#include "vtkTextProperty.h"
#include <vtkLookupTable.h>

DistanceMatrixLayer::DistanceMatrixLayer() {
	this->poly_data_ = vtkPolyData::New();
	this->mapper_ = vtkPolyDataMapper::New();
	this->actor_ = vtkActor::New();
	this->mapper_->SetInputData(this->poly_data_);
	this->actor_->SetMapper(this->mapper_);
	this->actor_->GetProperty()->SetLineWidth(2.0);
}

DistanceMatrixLayer::~DistanceMatrixLayer() {
}

void DistanceMatrixLayer::SetData(std::vector< std::string >& names, std::vector< std::vector< float > >& distance) {
	this->label_names_ = names;
	this->label_distance_ = distance;

	this->poly_data_->Initialize();

	vtkPoints* points = vtkPoints::New();
	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
	colors->SetNumberOfComponents(3);
	poly_data_->SetPoints(points);
	poly_data_->GetPointData()->SetScalars(colors);
	vtkCellArray* poly_array = vtkCellArray::New();
	poly_data_->SetPolys(poly_array);
	vtkCellArray* line_array = vtkCellArray::New();
	poly_data_->SetLines(line_array);

	float width_per_label = 0.7 / names.size();

	vtkIdType cell_ids[4];
	for (int i = 0; i < distance.size(); ++i)
		for (int j = 0; j < distance[i].size(); ++j) {
			float x = 0.20 + j * width_per_label;
			float y = 0.80 - i * width_per_label;
			cell_ids[0] = points->InsertNextPoint(x, y, 0);
			cell_ids[1] = points->InsertNextPoint(x + width_per_label, y, 0);
			cell_ids[2] = points->InsertNextPoint(x + width_per_label, y - width_per_label, 0);
			cell_ids[3] = points->InsertNextPoint(x, y - width_per_label, 0);

			if (distance[i][j] >= 0) {
				for (int k = 0; k < 4; ++k)
					colors->InsertNextTuple3(distance[i][j] * 255, distance[i][j] * 255, distance[i][j] * 255);
				poly_data_->InsertNextCell(VTK_POLYGON, 4, cell_ids);
			} else {
				for (int k = 0; k < 4; ++k)
					colors->InsertNextTuple3(128, 128, 128);

				vtkIdType temp_ids[2];
				temp_ids[0] = cell_ids[0];
				temp_ids[1] = cell_ids[2];
				poly_data_->InsertNextCell(VTK_LINE, 2, temp_ids);
				temp_ids[0] = cell_ids[1];
				temp_ids[1] = cell_ids[3];
				poly_data_->InsertNextCell(VTK_LINE, 2, temp_ids);
			}
			for (int k = 0; k < 4; ++k) {
				vtkIdType temp_ids[2];
				temp_ids[0] = cell_ids[k];
				temp_ids[1] = cell_ids[(k + 1) % 4];
				poly_data_->InsertNextCell(VTK_LINE, 2, temp_ids);
			}
		}
	this->actor_->Modified();

	title_actor_ = vtkTextActor::New();
	title_actor_->SetInput("Distance Matrix");

	if (hor_label_actors_.size() != names.size()) {
		int size = hor_label_actors_.size();
		for (int i = 0; i < names.size() - size; ++i) {
			vtkTextActor* temp_actor = vtkTextActor::New();
			hor_label_actors_.push_back(temp_actor);

			temp_actor = vtkTextActor::New();
			ver_label_actors_.push_back(temp_actor);
		}
	}

	for (int i = 0; i < names.size(); ++i) {
		hor_label_actors_[i]->SetInput(names[i].c_str());
		hor_label_actors_[i]->SetTextScaleModeToNone();
		hor_label_actors_[i]->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
		hor_label_actors_[i]->GetTextProperty()->SetFontSize(14);
		hor_label_actors_[i]->Modified();

		ver_label_actors_[i]->SetInput(names[i].c_str());
		ver_label_actors_[i]->SetTextScaleModeToNone();
		ver_label_actors_[i]->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
		ver_label_actors_[i]->GetTextProperty()->SetFontSize(14);
		ver_label_actors_[i]->Modified();
	}
}

void DistanceMatrixLayer::BuildRepresentation() {
	if (this->CurrentRenderer == NULL) return;
}

void DistanceMatrixLayer::SetEnabled(int enabling) {
	if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

	if (enabling) {
		if (this->Enabled) return;

		this->Enabled = 1;

		this->Interactor->AddObserver(vtkCommand::LeftButtonPressEvent,
			this->EventCallbackCommand, this->Priority);

		this->BuildRepresentation();

		float width_per_label = 0.7 / hor_label_actors_.size();
		double point_pos[4];
		for (int i = 0; i < hor_label_actors_.size(); ++i) {
			double* viewport = this->DefaultRenderer->GetViewport();
			DefaultRenderer->SetWorldPoint(0.2 + i * width_per_label, 0.01, 0, 1.0);
			DefaultRenderer->WorldToView();
			DefaultRenderer->GetViewPoint(point_pos);
			DefaultRenderer->ViewToNormalizedViewport(point_pos[0], point_pos[1], point_pos[2]);
			DefaultRenderer->NormalizedViewportToViewport(point_pos[0], point_pos[1]);
			hor_label_actors_[i]->SetPosition(point_pos[0], point_pos[1]);
			hor_label_actors_[i]->Modified();

			DefaultRenderer->SetWorldPoint(0.05, 0.1 + i * width_per_label, 0, 1.0);
			DefaultRenderer->WorldToView();
			DefaultRenderer->GetViewPoint(point_pos);
			DefaultRenderer->ViewToNormalizedViewport(point_pos[0], point_pos[1], point_pos[2]);
			DefaultRenderer->NormalizedViewportToViewport(point_pos[0], point_pos[1]);
			ver_label_actors_[i]->SetPosition(point_pos[0], point_pos[1]);
			ver_label_actors_[i]->Modified();
		}

		this->DefaultRenderer->AddActor(this->actor_);

		for (int i = 0; i < hor_label_actors_.size(); ++i) {
			this->DefaultRenderer->AddActor(hor_label_actors_[i]);
		}

		for (int i = 0; i < ver_label_actors_.size(); ++i) {
			this->DefaultRenderer->AddActor(ver_label_actors_[i]);
		}
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->Interactor->RemoveObserver(this->EventCallbackCommand);

		for (int i = 0; i < hor_label_actors_.size(); ++i) {
			this->DefaultRenderer->RemoveActor(hor_label_actors_[i]);
		}

		for (int i = 0; i < ver_label_actors_.size(); ++i) {
			this->DefaultRenderer->RemoveActor(ver_label_actors_[i]);
		}

		this->DefaultRenderer->RemoveActor(this->actor_);

		this->InvokeEvent(vtkCommand::DisableEvent, NULL);
	}

	this->Interactor->Render();
}

