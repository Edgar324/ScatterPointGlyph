#include "transmap.h"
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
#include "transmap_data.h"
#include "vtkCellArray.h"

vtkStandardNewMacro(TransMap);

TransMap::TransMap() {
	this->state = TransMap::START;
	this->EventCallbackCommand->SetCallback(TransMap::ProcessEvents);

	float picker_tolerence = 1e-10;
	this->node_picker = vtkCellPicker::New();
	this->node_picker->SetTolerance(picker_tolerence);
	
	this->trans_picker = vtkCellPicker::New();
	this->trans_picker->SetTolerance(picker_tolerence);

	this->boundary_picker = vtkCellPicker::New();
	this->boundary_picker->SetTolerance(picker_tolerence);

	this->current_handle = NULL;

	double bounds[6];
	bounds[0] = -0.5;
	bounds[1] = 0.5;
	bounds[2] = -0.5;
	bounds[3] = 0.5;
	bounds[4] = -0.5;
	bounds[5] = 0.5;
	this->PlaceFactor = 1.0;
	this->PlaceWidget(bounds);

	this->node_radius = 10;
}

TransMap::~TransMap() {
	
}

void TransMap::SetData(TransMapData* data) {
	dataset_ = data;

	this->ClearActors();
	this->InitActors();
	this->BuildRepresentation();
}

void TransMap::InitActors() {

	vtkIdType cell_ids[5];

	for (int i = 0; i < dataset_->node_num; ++i) {
		// construct the polydata
		vtkPoints* points = vtkPoints::New();
		vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
		colors->SetNumberOfComponents(3);

		vtkPolyData* polydata = vtkPolyData::New();
		polydata->SetPoints(points);
		polydata->GetPointData()->SetScalars(colors);

		vtkCellArray* poly_array = vtkCellArray::New();
		polydata->SetPolys(poly_array);

		for (int j = 0; j < dataset_->var_num; ++j) {
			float temp_radius = node_radius * dataset_->node_average_value[i][j] * 0.8 + node_radius * 0.2;
			vtkIdType center_id = points->InsertNextPoint(dataset_->node_center[i][0], dataset_->node_center[i][1], 0);

			float end_arc = j * 3.14159 * 2 / dataset_->var_num;
			float x = temp_radius * cos(end_arc);
			float y = temp_radius * sin(end_arc);

			vtkIdType pre_id = points->InsertNextPoint(dataset_->node_center[i][0] + x, dataset_->node_center[i][1] + y, 0);
			for (int k = 1; k <= 8; ++k) {
				end_arc = ((float)k / 8 + j) * 3.14159 * 2 / dataset_->var_num;
				x = temp_radius * cos(end_arc);
				y = temp_radius * sin(end_arc);
				vtkIdType current_id = points->InsertNextPoint(dataset_->node_center[i][0] + x, dataset_->node_center[i][1] + y, 0);

				cell_ids[0] = center_id;
				cell_ids[1] = pre_id;
				cell_ids[2] = current_id;
				polydata->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

				pre_id = current_id;
			}

			int rgb[3];
			rgb[0] = dataset_->var_repsentative_color[3 * i];
			rgb[1] = dataset_->var_repsentative_color[3 * i + 1];
			rgb[2] = dataset_->var_repsentative_color[3 * i + 2];
			for (int k = 0; k < 10; ++k) colors->InsertNextTuple3(rgb[0], rgb[1], rgb[2]);
		}

		vtkActor* node_actor = vtkActor::New();
		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInputData(polydata);
		node_actor->SetMapper(mapper);
		
		node_glyph_actors.push_back(node_actor);
	}

	for (int i = 0; i < dataset_->node_num - 1; ++i)
		for (int j = i + 1; j < dataset_->node_num; ++j) 
			if (dataset_->node_connecting_status[i][j]) {

				vtkPoints* points = vtkPoints::New();
				vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
				colors->SetNumberOfComponents(3);

				vtkPolyData* polydata = vtkPolyData::New();
				polydata->SetPoints(points);
				polydata->GetPointData()->SetScalars(colors);

				vtkCellArray* poly_array = vtkCellArray::New();
				polydata->SetPolys(poly_array);

				float band_width = 1.5 * node_radius;
				float width_per_var = band_width / dataset_->var_num;
				float dir_vec[2], orth_vec[2];
				dir_vec[0] = dataset_->node_center[j][0] - dataset_->node_center[i][0];
				dir_vec[1] = dataset_->node_center[j][1] - dataset_->node_center[i][1];
				float dir_length = sqrt(pow(dir_vec[0], 2) + pow(dir_vec[1], 2));
				dir_vec[0] /= dir_length;
				dir_vec[1] /= dir_length;
				orth_vec[0] = -1 * dir_vec[1];
				orth_vec[1] = dir_vec[0];

				float begin_pos[2], end_pos[2];
				begin_pos[0] = dataset_->node_center[i][0] + node_radius * dir_vec[0];
				begin_pos[1] = dataset_->node_center[i][1] + node_radius * dir_vec[1];
				end_pos[0] = dataset_->node_center[j][0] - node_radius * dir_vec[0];
				end_pos[1] = dataset_->node_center[j][1] - node_radius * dir_vec[1];

				int var_index = 2;

				cell_ids[0] = points->InsertNextPoint(begin_pos[0] + 0.5 * node_radius * orth_vec[0] * dataset_->node_average_value[i][var_index], begin_pos[1] + 0.5 * node_radius * orth_vec[1] * dataset_->node_average_value[i][var_index], 0);
				cell_ids[1] = points->InsertNextPoint(end_pos[0] + 0.5 * node_radius * orth_vec[0] * dataset_->node_average_value[j][var_index], end_pos[1] + 0.5 * node_radius * orth_vec[1] * dataset_->node_average_value[j][var_index], 0);
				cell_ids[2] = points->InsertNextPoint(end_pos[0], end_pos[1], 0);
				cell_ids[3] = points->InsertNextPoint(begin_pos[0], begin_pos[1], 0);

				polydata->InsertNextCell(VTK_POLYGON, 4, cell_ids);

				/*cell_ids[0] = points->InsertNextPoint(begin_pos[0], begin_pos[1], 0);
				cell_ids[1] = points->InsertNextPoint(end_pos[0], end_pos[1], 0);
				cell_ids[2] = points->InsertNextPoint(end_pos[0] - 0.5 * node_radius * orth_vec[0] * dataset_->node_average_value[j][var_index + 1], end_pos[1] - 0.5 * node_radius * orth_vec[1] * dataset_->node_average_value[j][var_index + 1], 0);
				cell_ids[3] = points->InsertNextPoint(begin_pos[0] - 0.5 * node_radius * orth_vec[0] * dataset_->node_average_value[i][var_index + 1], begin_pos[1] - 0.5 * node_radius * orth_vec[1] * dataset_->node_average_value[i][var_index + 1], 0);
				polydata->InsertNextCell(VTK_POLYGON, 4, cell_ids);*/

				int color = 200 * dataset_->node_average_value[i][var_index];
				int color2 = 200 * dataset_->node_average_value[j][var_index];
				colors->InsertNextTuple3(color, 0, 0);
				colors->InsertNextTuple3(color2, 0, 0);
				colors->InsertNextTuple3(color2, 0, 0);
				colors->InsertNextTuple3(color, 0, 0);

				/*color = 200 * dataset_->node_average_value[i][var_index + 1];
				color2 = 200 * dataset_->node_average_value[j][var_index + 1];
				colors->InsertNextTuple3(0, color, 0);
				colors->InsertNextTuple3(0, color2, 0);
				colors->InsertNextTuple3(0, color2, 0);
				colors->InsertNextTuple3(0, color, 0);*/

				vtkActor* node_actor = vtkActor::New();
				vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
				mapper->SetInputData(polydata);
				node_actor->SetMapper(mapper);

				trans_glyph_actors.push_back(node_actor);
			}
}

void TransMap::ClearActors() {
	for (int i = 0; i < node_glyph_actors.size(); ++i) {
		this->CurrentRenderer->RemoveActor(node_glyph_actors[i]);
	}
	node_glyph_actors.clear();

	for (int i = 0; i < trans_glyph_actors.size(); ++i) {
		this->CurrentRenderer->RemoveActor(trans_glyph_actors[i]);
	}
	trans_glyph_actors.clear();
}

void TransMap::SetEnabled(int enabling) {
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
	}

	vtkRenderWindowInteractor *i = this->Interactor;
	i->AddObserver(vtkCommand::MouseMoveEvent,
		this->EventCallbackCommand, this->Priority);
	i->AddObserver(vtkCommand::LeftButtonPressEvent,
		this->EventCallbackCommand, this->Priority);
	i->AddObserver(vtkCommand::LeftButtonReleaseEvent,
		this->EventCallbackCommand, this->Priority);
	i->AddObserver(vtkCommand::RightButtonPressEvent,
		this->EventCallbackCommand, this->Priority);
	i->AddObserver(vtkCommand::RightButtonReleaseEvent,
		this->EventCallbackCommand, this->Priority);

	for (int i = 0; i < node_glyph_actors.size(); ++i) {
		this->CurrentRenderer->AddActor(node_glyph_actors[i]);
	}

	for (int i = 0; i < trans_glyph_actors.size(); ++i) {
		this->CurrentRenderer->AddActor(trans_glyph_actors[i]);
	}

	this->Interactor->Render();
}

void TransMap::BuildRepresentation() {
	
}

void TransMap::ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata) {
	TransMap* self = reinterpret_cast< TransMap* >(clientdata);

	switch (event)
	{
	case vtkCommand::LeftButtonPressEvent:
		self->OnLeftButtonDown();
		break;
	case vtkCommand::LeftButtonReleaseEvent:
		self->OnLeftButtonUp();
		break;
	case vtkCommand::RightButtonPressEvent:
		self->OnRightButtonDown();
		break;
	case vtkCommand::RightButtonReleaseEvent:
		self->OnRightButtonUp();
		break;
	case vtkCommand::MouseMoveEvent:
		self->OnMouseMove();
		break;
	}
}

void TransMap::OnLeftButtonDown() {
}

void TransMap::OnLeftButtonUp() {
}

void TransMap::OnRightButtonDown() {
}

void TransMap::OnRightButtonUp() {
}

void TransMap::OnMouseMove() {
}
