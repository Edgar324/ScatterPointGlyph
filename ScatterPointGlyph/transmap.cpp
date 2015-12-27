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
#include <QVTKWidget.h>
#include "transmap_data.h"
#include "vtkCellArray.h"
#include "tree_common.h"

TransMap::TransMap() {
	this->state = TransMap::NORMAL;
	this->EventCallbackCommand->SetCallback(TransMap::ProcessEvents);

	float picker_tolerence = 1e-10;
	this->level_one_node_picker = vtkCellPicker::New();
	this->level_one_node_picker->SetTolerance(picker_tolerence);
	this->level_one_node_picker->PickFromListOn();
	
	this->level_zero_node_picker = vtkCellPicker::New();
	this->level_zero_node_picker->SetTolerance(picker_tolerence);
	this->level_zero_node_picker->PickFromListOn();

	this->trans_picker = vtkCellPicker::New();
	this->trans_picker->SetTolerance(picker_tolerence);
	this->trans_picker->PickFromListOn();

	this->boundary_picker = vtkCellPicker::New();
	this->boundary_picker->SetTolerance(picker_tolerence);
	this->boundary_picker->PickFromListOn();

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
	this->is_sequence_selection = true;
	this->current_handle = NULL;
	this->current_node = NULL;

	this->highlight_actor = vtkActor::New();
	this->highlight_poly = vtkPolyData::New();
	this->hightlight_mapper = vtkPolyDataMapper::New();
	this->hightlight_mapper->SetInputData(this->highlight_poly);
	this->highlight_actor->SetMapper(this->hightlight_mapper);
	this->highlight_actor->GetProperty()->SetLineWidth(3.0);
	

	this->selection_brush_poly = vtkPolyData::New();
	this->selection_brush_mapper = vtkPolyDataMapper::New();
	this->selection_brush_actor = vtkActor::New();
	this->selection_brush_mapper->SetInputData(this->selection_brush_poly);
	this->selection_brush_actor->SetMapper(this->selection_brush_mapper);
	this->selection_brush_actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
	this->selection_brush_actor->GetProperty()->SetLineWidth(3.0);
}

TransMap::TransMap(QVTKWidget* parent) 
	: TransMap() {
	this->parent_view = parent;
}

TransMap::~TransMap() {
	
}

void TransMap::SetData(TransMapData* data) {
	dataset_ = data;

	this->SetEnabled(false);
	this->ClearActors();
	this->ConstructActors();
	this->BuildRepresentation();
	this->SetEnabled(true);

	highlight_node_sequence.clear();
}

void TransMap::SetSequenceSelectionOn() {
	this->is_sequence_selection = true;
}

void TransMap::SetSequenceSelectionOff()  {
	this->is_sequence_selection = false;
}

void TransMap::SetBrushSelectionOn() {
	this->state_vec.push_back(this->state);
	this->state = WidgetState::SELECT_CLUSTER;
	this->highlight_actor->SetVisibility(false);
}

void TransMap::SetBrushSelectionOff() {
	this->state = this->state_vec[this->state_vec.size() - 1];
	this->state_vec.pop_back();
	this->selection_brush_poly->Initialize();
	this->selection_brush_actor->Modified();
	this->highlight_actor->SetVisibility(true);
}

void TransMap::SetMouseReleased() {
	this->OnLeftButtonUp();
}

void TransMap::SetMouseDragmove(int x, int y) {
	if (this->state == WidgetState::SELECTION_BRUSH_BEGIN) {
		double display_pos[4];
		display_pos[0] = x;
		display_pos[1] = y;
		display_pos[2] = 0;
		display_pos[3] = 1;
		this->DefaultRenderer->SetDisplayPoint(display_pos);
		this->DefaultRenderer->DisplayToWorld();
		double* world_pos = this->DefaultRenderer->GetWorldPoint();

		this->selection_brush_poly->GetPoints()->InsertNextPoint(world_pos[0], world_pos[1], 0);

		vtkCellArray* line_array = this->selection_brush_poly->GetLines();
		line_array->Initialize();
		int pnum = this->selection_brush_poly->GetPoints()->GetNumberOfPoints();
		vtkIdType cell_ids[2];
		for (int i = 0; i < pnum; ++i){
			cell_ids[0] = i;
			cell_ids[1] = (i + 1) % pnum;
			this->selection_brush_poly->InsertNextCell(VTK_LINE, 2, cell_ids);
		}

		this->selection_brush_poly->Modified();
		this->selection_brush_actor->Modified();

		this->parent_view->update();
	}
}

int TransMap::GetSelectedClusterIndex() {
	if (this->current_node != NULL) {
		std::map< int, CNode* >::iterator iter = dataset_->cluster_node_map.begin();
		while (iter != dataset_->cluster_node_map.end() && iter->second != this->current_node) iter++;
		if (iter != dataset_->cluster_node_map.end()) {
			return iter->first;
		}
	} 
	return -1;
}

void TransMap::GetSelectedClusterIndex(std::vector< int >& index) {
	index.clear();

	std::list< CNode* >::iterator node_iter = highlight_node_sequence.begin();
	while (node_iter != highlight_node_sequence.end()) {
		std::map< int, CNode* >::iterator iter = dataset_->cluster_node_map.begin();
		while (iter != dataset_->cluster_node_map.end() && iter->second != *node_iter) iter++;
		if (iter != dataset_->cluster_node_map.end()) {
			index.push_back(iter->first);
		}
		node_iter++;
	}
}

void TransMap::ConstructActors() {
	vtkIdType cell_ids[5];

	for (int i = 0; i < dataset_->level_one_nodes.size(); ++i) {
		CNode* node = dataset_->level_one_nodes[i];
		// construct the polydata
		vtkPoints* points = vtkPoints::New();
		vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
		colors->SetNumberOfComponents(4);

		vtkPolyData* polydata = vtkPolyData::New();
		polydata->SetPoints(points);
		polydata->GetPointData()->SetScalars(colors);

		vtkCellArray* poly_array = vtkCellArray::New();
		polydata->SetPolys(poly_array);

		vtkCellArray* line_array = vtkCellArray::New();
		polydata->SetLines(line_array);

		// insert background
		std::vector< vtkIdType > background_ids;
		for (int j = 0; j < dataset_->var_num; ++j) {
			float end_arc = j * 3.14159 * 2 / dataset_->var_num;
			float x = node_radius * cos(end_arc);
			float y = node_radius * sin(end_arc);

			background_ids.push_back(points->InsertNextPoint(node->center_pos[0] + x, node->center_pos[1] + y, -0.00002));
			colors->InsertNextTuple4(128, 128, 128, 10);
		}
		polydata->InsertNextCell(VTK_POLYGON, dataset_->var_num, background_ids.data());

		// paint bound
		std::vector< int > var_point_ids;
		var_point_ids.resize(dataset_->var_num * 2);
		for (int j = 0; j < dataset_->var_num; ++j) {
			float end_arc = j * 3.14159 * 2 / dataset_->var_num;

			float temp_radius = node_radius * node->boxplot_lower_bound[j]/* * 0.8 + node_radius * 0.2*/;
			float x = temp_radius * cos(end_arc);
			float y = temp_radius * sin(end_arc);

			vtkIdType id_one = points->InsertNextPoint(node->center_pos[0] + x, node->center_pos[1] + y, -0.00001);
			colors->InsertNextTuple4(128, 128, 128, 128);

			temp_radius = node_radius * node->boxplot_upper_bound[j]/* * 0.8 + node_radius * 0.2*/;
			x = temp_radius * cos(end_arc);
			y = temp_radius * sin(end_arc);

			vtkIdType id_two = points->InsertNextPoint(node->center_pos[0] + x, node->center_pos[1] + y, -0.00001);
			colors->InsertNextTuple4(128, 128, 128, 128);

			var_point_ids[j * 2] = id_one;
			var_point_ids[j * 2 + 1] = id_two;
		}
		for (int j = 0; j < dataset_->var_num; ++j) {
			cell_ids[0] = var_point_ids[2 * j];
			cell_ids[1] = var_point_ids[2 * j + 1];
			cell_ids[2] = var_point_ids[2 * ((j + 1) % dataset_->var_num)];
			polydata->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

			cell_ids[0] = var_point_ids[2 * j + 1];
			cell_ids[1] = var_point_ids[2 * ((j + 1) % dataset_->var_num) + 1];
			cell_ids[2] = var_point_ids[2 * ((j + 1) % dataset_->var_num)];
			polydata->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
		}

		// paint value
		std::vector< vtkIdType > value_pos_ids;
		value_pos_ids.resize(dataset_->var_num);
		for (int j = 0; j < dataset_->var_num; ++j) {
			float temp_radius = node_radius * node->average_values[j]/* * 0.8 + node_radius * 0.2*/;
			float end_arc = j * 3.14159 * 2 / dataset_->var_num;
			float x = temp_radius * cos(end_arc);
			float y = temp_radius * sin(end_arc);

			value_pos_ids[j] = points->InsertNextPoint(node->center_pos[0] + x, node->center_pos[1] + y, 0);
			colors->InsertNextTuple4(dataset_->level_one_colors[3 * i], dataset_->level_one_colors[3 * i + 1], dataset_->level_one_colors[3 * i + 2], 255);
		}
		for (int j = 0; j < dataset_->var_num; ++j) {
			cell_ids[0] = value_pos_ids[j];
			cell_ids[1] = value_pos_ids[(j + 1) % dataset_->var_num];
			polydata->InsertNextCell(VTK_LINE, 2, cell_ids);
		}
		

		// paint axis
		vtkIdType center_id = points->InsertNextPoint(node->center_pos[0], node->center_pos[1], -0.0001);
		colors->InsertNextTuple4(0, 0, 0, 255);
		for (int j = 0; j < dataset_->var_num; ++j) {
			float end_arc = j * 3.14159 * 2 / dataset_->var_num;
			float x = node_radius * cos(end_arc);
			float y = node_radius * sin(end_arc);

			vtkIdType pre_id = points->InsertNextPoint(node->center_pos[0] + x, node->center_pos[1] + y, -0.0001);
			colors->InsertNextTuple4(0, 0, 0, 255);
			
			cell_ids[0] = center_id;
			cell_ids[1] = pre_id;
			polydata->InsertNextCell(VTK_LINE, 2, cell_ids);
		}

		vtkActor* node_actor = vtkActor::New();
		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInputData(polydata);
		node_actor->SetMapper(mapper);
		node_actor->GetProperty()->SetLineWidth(2.0);
		
		level_one_node_glyph_actors.push_back(node_actor);
	}

	for (int i = 0; i < dataset_->level_zero_nodes.size(); ++i) {
		CNode* node = dataset_->level_zero_nodes[i];
		// construct the polydata
		vtkPoints* points = vtkPoints::New();
		vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
		colors->SetNumberOfComponents(3);

		vtkPolyData* polydata = vtkPolyData::New();
		polydata->SetPoints(points);
		polydata->GetPointData()->SetScalars(colors);

		vtkCellArray* poly_array = vtkCellArray::New();
		polydata->SetPolys(poly_array);

		float temp_radius = node_radius * 0.2;
		float end_arc = 0;
		float x = temp_radius;
		float y = 0;

		vtkIdType center_id = points->InsertNextPoint(node->center_pos[0], node->center_pos[1], 0);
		vtkIdType pre_id = points->InsertNextPoint(node->center_pos[0] + x, node->center_pos[1] + y, 0);

		for (int k = 1; k <= 20; ++k) {
			end_arc = (float)k / 20 * 3.14159 * 2;
			x = temp_radius * cos(end_arc);
			y = temp_radius * sin(end_arc);
			vtkIdType current_id = points->InsertNextPoint(node->center_pos[0] + x, node->center_pos[1] + y, 0);

			cell_ids[0] = pre_id;
			cell_ids[1] = current_id;
			cell_ids[2] = center_id;
			polydata->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

			pre_id = current_id;
		}

		for (int k = 0; k < 22; ++k) colors->InsertNextTuple3(dataset_->level_zero_colors[3 * i], dataset_->level_zero_colors[3 * i + 1], dataset_->level_zero_colors[3 * i + 2]);

		vtkActor* node_actor = vtkActor::New();
		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInputData(polydata);
		node_actor->SetMapper(mapper);
		node_actor->GetProperty()->SetLineWidth(2.0);

		level_zero_node_glyph_actors.push_back(node_actor);
	}

	for (int i = 0; i < dataset_->level_one_nodes.size() - 1; ++i){
		CNode* node_one = dataset_->level_one_nodes[i];

		for (int j = i + 1; j < dataset_->level_one_nodes.size(); ++j) {
			if (dataset_->node_connecting_status[i][j]) {
				CNode* node_two = dataset_->level_one_nodes[j];

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
				dir_vec[0] = node_two->center_pos[0] - node_one->center_pos[0];
				dir_vec[1] = node_two->center_pos[1] - node_one->center_pos[1];
				float dir_length = sqrt(pow(dir_vec[0], 2) + pow(dir_vec[1], 2));
				dir_vec[0] /= dir_length;
				dir_vec[1] /= dir_length;
				orth_vec[0] = -1 * dir_vec[1];
				orth_vec[1] = dir_vec[0];

				float begin_pos[2], end_pos[2];
				begin_pos[0] = node_one->center_pos[0] + node_radius * dir_vec[0];
				begin_pos[1] = node_one->center_pos[1] + node_radius * dir_vec[1];
				end_pos[0] = node_two->center_pos[0] - node_radius * dir_vec[0];
				end_pos[1] = node_two->center_pos[1] - node_radius * dir_vec[1];
				
				float value_dis = 0;
				for (int k = 0; k < dataset_->var_num; ++k)
					value_dis += abs(node_one->average_values[k] - node_two->average_values[k]);
				value_dis /= dataset_->var_num;

				cell_ids[0] = points->InsertNextPoint(begin_pos[0] + 0.5 * node_radius * orth_vec[0] * value_dis, begin_pos[1] + 0.5 * node_radius * orth_vec[1] * value_dis, 0);
				cell_ids[1] = points->InsertNextPoint(end_pos[0] + 0.5 * node_radius * orth_vec[0] * value_dis, end_pos[1] + 0.5 * node_radius * orth_vec[1] * value_dis, 0);
				cell_ids[2] = points->InsertNextPoint(end_pos[0], end_pos[1], 0);
				cell_ids[3] = points->InsertNextPoint(begin_pos[0], begin_pos[1], 0);

				polydata->InsertNextCell(VTK_POLYGON, 4, cell_ids);

				colors->InsertNextTuple3(128, 128, 128);
				colors->InsertNextTuple3(128, 128, 128);
				colors->InsertNextTuple3(128, 128, 128);
				colors->InsertNextTuple3(128, 128, 128);

				vtkActor* node_actor = vtkActor::New();
				vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
				mapper->SetInputData(polydata);
				node_actor->SetMapper(mapper);

				trans_glyph_actors.push_back(node_actor);
			}
		}
	}
}

void TransMap::ClearActors() {
	for (int i = 0; i < level_one_node_glyph_actors.size(); ++i) {
		this->level_one_node_picker->DeletePickList(level_one_node_glyph_actors[i]);
		level_one_node_glyph_actors[i]->Delete();
	}
	level_one_node_glyph_actors.clear();

	for (int i = 0; i < level_zero_node_glyph_actors.size(); ++i) {
		this->level_zero_node_picker->DeletePickList(level_zero_node_glyph_actors[i]);
		level_zero_node_glyph_actors[i]->Delete();
	}
	level_zero_node_glyph_actors.clear();

	for (int i = 0; i < trans_glyph_actors.size(); ++i) {
		this->trans_picker->DeletePickList(trans_glyph_actors[i]);
		trans_glyph_actors[i]->Delete();
	}
	trans_glyph_actors.clear();

	this->highlight_poly->Initialize();
	this->highlight_actor->Modified();
}

void TransMap::SetEnabled(int enabling) {
	if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

	if (this->DefaultRenderer) this->SetCurrentRenderer(this->DefaultRenderer);

	if (enabling) {
		if (this->Enabled) return;

		this->Enabled = 1;

		this->Interactor->AddObserver(vtkCommand::LeftButtonPressEvent,
			this->EventCallbackCommand, this->Priority);
		this->Interactor->AddObserver(vtkCommand::MouseMoveEvent,
			this->EventCallbackCommand, this->Priority);
		this->Interactor->AddObserver(vtkCommand::LeftButtonReleaseEvent,
			this->EventCallbackCommand, this->Priority);
		this->Interactor->AddObserver(vtkCommand::RightButtonPressEvent,
			this->EventCallbackCommand, this->Priority);
		this->Interactor->AddObserver(vtkCommand::RightButtonReleaseEvent,
			this->EventCallbackCommand, this->Priority);

		for (int i = 0; i < level_one_node_glyph_actors.size(); ++i) {
			this->DefaultRenderer->AddActor(level_one_node_glyph_actors[i]);
			this->level_one_node_picker->AddPickList(level_one_node_glyph_actors[i]);
		}

		for (int i = 0; i < level_zero_node_glyph_actors.size(); ++i) {
			this->DefaultRenderer->AddActor(level_zero_node_glyph_actors[i]);
			this->level_zero_node_picker->AddPickList(level_zero_node_glyph_actors[i]);
		}

		for (int i = 0; i < trans_glyph_actors.size(); ++i) {
			this->DefaultRenderer->AddActor(trans_glyph_actors[i]);
			this->trans_picker->AddPickList(trans_glyph_actors[i]);
		}

		this->DefaultRenderer->AddActor(this->highlight_actor);
		this->DefaultRenderer->AddActor(this->selection_brush_actor);

		this->InvokeEvent(vtkCommand::EnableEvent, NULL);
	} else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->Interactor->RemoveObserver(this->EventCallbackCommand);

		for (int i = 0; i < level_one_node_glyph_actors.size(); ++i) {
			this->DefaultRenderer->RemoveActor(level_one_node_glyph_actors[i]);
			this->level_one_node_picker->DeletePickList(level_one_node_glyph_actors[i]);
		}

		for (int i = 0; i < level_zero_node_glyph_actors.size(); ++i) {
			this->DefaultRenderer->RemoveActor(level_zero_node_glyph_actors[i]);
			this->level_zero_node_picker->DeletePickList(level_zero_node_glyph_actors[i]);
		}

		for (int i = 0; i < trans_glyph_actors.size(); ++i) {
			this->DefaultRenderer->RemoveActor(trans_glyph_actors[i]);
			this->trans_picker->DeletePickList(trans_glyph_actors[i]);
		}

		this->DefaultRenderer->RemoveActor(this->highlight_actor);
		this->DefaultRenderer->RemoveActor(this->selection_brush_actor);

		this->Interactor->RemoveObserver(vtkCommand::LeftButtonPressEvent);
		this->Interactor->RemoveObserver(vtkCommand::MouseMoveEvent);
		this->Interactor->RemoveObserver(vtkCommand::LeftButtonReleaseEvent);
		this->Interactor->RemoveObserver(vtkCommand::RightButtonPressEvent);
		this->Interactor->RemoveObserver(vtkCommand::RightButtonReleaseEvent);

		this->InvokeEvent(vtkCommand::DisableEvent, NULL);
		this->SetCurrentRenderer(NULL);
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
	default:
		break;
	}
}

void TransMap::OnLeftButtonDown() {
	int X = this->Interactor->GetEventPosition()[0];
	int Y = this->Interactor->GetEventPosition()[1];

	if (!this->CurrentRenderer) {
		this->state = TransMap::NORMAL;
		this->current_handle = NULL;
		this->current_node = NULL;
		this->highlight_poly->Initialize();
		this->highlight_actor->Modified();
		return;
	}

	if (this->state == WidgetState::SELECT_CLUSTER) {
		this->state_vec.push_back(this->state);
		this->state = WidgetState::SELECTION_BRUSH_BEGIN;

		this->selection_brush_poly->Initialize();

		vtkPoints* points = vtkPoints::New();
		selection_brush_poly->SetPoints(points);

		vtkCellArray* line_array = vtkCellArray::New();
		selection_brush_poly->SetLines(line_array);

		this->selection_brush_actor->Modified();

		this->parent_view->update();

		return;
	}

	vtkAssemblyPath* path = this->GetAssemblyPath(X, Y, 0., this->level_one_node_picker);
	if (path != NULL) {
		this->state = TransMap::HIGHLIGHT_LEVEL_ONE_NODE;
		this->HighlightHandle(path->GetFirstNode()->GetViewProp());
		return;
	}

	path = this->GetAssemblyPath(X, Y, 0., this->level_zero_node_picker);
	if (path != NULL) {
		this->state = TransMap::HIGHLIGHT_LEVEL_ZERO_NODE;
		this->HighlightHandle(path->GetFirstNode()->GetViewProp());
		return;
	}

	path = this->GetAssemblyPath(X, Y, 0., this->trans_picker);
	if (path != NULL) {
		this->state = TransMap::HIGHLIGHT_TRANSFER_BAND;
		this->HighlightHandle(path->GetFirstNode()->GetViewProp());
		return;
	}
}

void TransMap::OnLeftButtonUp() {
	if (this->state == WidgetState::SELECTION_BRUSH_BEGIN) {
		this->state = this->state_vec[this->state_vec.size() - 1];
		this->state_vec.pop_back();

		this->UpdateBrushSelectionResult();

		this->selection_brush_poly->Initialize();
		this->selection_brush_actor->Modified();
	}
}

void TransMap::OnRightButtonDown() {
}

void TransMap::OnRightButtonUp() {
}

void TransMap::OnMouseMove() {
	if (this->state == WidgetState::SELECTION_BRUSH_BEGIN) {
		
	}
}

void TransMap::HighlightHandle(vtkProp* prop) {
	if (prop == NULL) return;

	vtkActor* actor = dynamic_cast< vtkActor* >(prop);
	//if (this->current_handle == actor) return;
	if (actor == NULL) return;

	switch (this->state) {
	case TransMap::HIGHLIGHT_LEVEL_ZERO_NODE:
	{
		int node_index = -1;
		for (int i = 0; i < level_zero_node_glyph_actors.size(); ++i)
			if (level_zero_node_glyph_actors[i] == actor) {
				node_index = i;
				break;
			}
		if (node_index == -1) break;
		this->current_node = dataset_->level_zero_nodes[node_index];
		this->SelectNode(this->current_node);
	}
		break;
	case TransMap::HIGHLIGHT_LEVEL_ONE_NODE:
	{
		int node_index = -1;
		for (int i = 0; i < level_one_node_glyph_actors.size(); ++i)
			if (level_one_node_glyph_actors[i] == actor) {
				node_index = i;
				break;
			}
		if (node_index == -1) break;
		this->current_node = dataset_->level_one_nodes[node_index];
		this->SelectNode(this->current_node);
	}
		break;
	default:
		break;
	}

	this->UpdateHightlightActor();
}

void TransMap::UpdateHightlightActor() {
	this->highlight_poly->Initialize();

	vtkPoints* points = vtkPoints::New();
	vtkUnsignedCharArray* color_array = vtkUnsignedCharArray::New();
	color_array->SetNumberOfComponents(3);

	highlight_poly->SetPoints(points);
	highlight_poly->GetPointData()->SetScalars(color_array);

	vtkCellArray* line_array = vtkCellArray::New();
	highlight_poly->SetLines(line_array);

	float end_arc = 0;
	float x = node_radius;
	float y = 0;
	vtkIdType cell_ids[5];

	std::list< CNode* >::iterator iter = highlight_node_sequence.begin();
	while (iter != highlight_node_sequence.end()) {
		int node_index = -1;
		for (int i = 0; i < dataset_->level_one_nodes.size(); ++i)
			if (dataset_->level_one_nodes[i] == *iter) {
				node_index = i;
				break;
			}
		if (node_index != -1) {
			CNode* temp_node = dataset_->level_one_nodes[node_index];
			vtkIdType pre_id = points->InsertNextPoint(temp_node->center_pos[0] + x, temp_node->center_pos[1] + y, 0);
			for (int k = 1; k <= 20; ++k) {
				end_arc = (float)k / 20 * 3.14159 * 2;
				x = node_radius * cos(end_arc);
				y = node_radius * sin(end_arc);
				vtkIdType current_id = points->InsertNextPoint(temp_node->center_pos[0] + x, temp_node->center_pos[1] + y, 0);

				cell_ids[0] = pre_id;
				cell_ids[1] = current_id;
				this->highlight_poly->InsertNextCell(VTK_LINE, 2, cell_ids);

				pre_id = current_id;
			}

			for (int k = 0; k < 21; ++k) color_array->InsertNextTuple3(255, 0, 0);
			iter++;
			continue;
		}

		for (int i = 0; i < dataset_->level_zero_nodes.size(); ++i)
			if (dataset_->level_zero_nodes[i] == *iter) {
				node_index = i;
				break;
			}
		if (node_index != -1) {
			CNode* temp_node = dataset_->level_zero_nodes[node_index];
			vtkIdType pre_id = points->InsertNextPoint(temp_node->center_pos[0] + x * 0.4, temp_node->center_pos[1] + y * 0.4, 0);
			for (int k = 1; k <= 20; ++k) {
				end_arc = (float)k / 20 * 3.14159 * 2;
				x = node_radius * cos(end_arc);
				y = node_radius * sin(end_arc);
				vtkIdType current_id = points->InsertNextPoint(temp_node->center_pos[0] + x * 0.4, temp_node->center_pos[1] + y * 0.4, 0);

				cell_ids[0] = pre_id;
				cell_ids[1] = current_id;
				this->highlight_poly->InsertNextCell(VTK_LINE, 2, cell_ids);

				pre_id = current_id;
			}

			for (int k = 0; k < 21; ++k) color_array->InsertNextTuple3(255, 0, 0);
			iter++;
			continue;
		}

		iter++;
	}

	this->highlight_actor->GetProperty()->SetLineStipplePattern(0xF00F);
	this->highlight_actor->Modified();
}

void TransMap::SelectNode(CNode* node) {
	if (!is_sequence_selection) {
		if (highlight_node_sequence.size() > 1) {
			highlight_node_sequence.clear();
			highlight_node_sequence.push_back(node);
		} else if (highlight_node_sequence.size() == 1) {
			if (*highlight_node_sequence.begin() == node) {
				highlight_node_sequence.clear();
			}
			else {
				highlight_node_sequence.clear();
				highlight_node_sequence.push_back(node);
			}
		} else {
			highlight_node_sequence.push_back(node);
		}
	} else {
		std::list< CNode* >::iterator iter = highlight_node_sequence.begin();
		while (iter != highlight_node_sequence.end() && *iter != node) iter++;
		if (iter != highlight_node_sequence.end()) {
			highlight_node_sequence.erase(iter);
		} else {
			highlight_node_sequence.push_back(node);
		}
	}
}

void TransMap::UpdateBrushSelectionResult() {
	highlight_node_sequence.clear();

	for (int i = 0; i < dataset_->level_one_nodes.size(); ++i) {
		if (IsInsideSelection(dataset_->level_one_nodes[i]->center_pos[0], dataset_->level_one_nodes[i]->center_pos[1]))
			highlight_node_sequence.push_back(dataset_->level_one_nodes[i]);
	}

	for (int i = 0; i < dataset_->level_zero_nodes.size(); ++i) {
		if (IsInsideSelection(dataset_->level_zero_nodes[i]->center_pos[0], dataset_->level_zero_nodes[i]->center_pos[1]))
			highlight_node_sequence.push_back(dataset_->level_zero_nodes[i]);
	}

	this->UpdateHightlightActor();
}

bool TransMap::IsInsideSelection(float x, float y) {
	int count = this->selection_brush_poly->GetPoints()->GetNumberOfPoints();

	if (count < 3) return false;

	const double PI = 3.14159265358979323846;
	double angle = 0;

	for (int i = 0, k = count - 1; i < count; k = i++) {
		double p1[2];
		double p2[2];
		double* temp_pos = this->selection_brush_poly->GetPoint(i);
		p1[0] = temp_pos[0] - x;
		p1[1] = temp_pos[1] - y;

		temp_pos = this->selection_brush_poly->GetPoint(k);
		p2[0] = temp_pos[0] - x;
		p2[1] = temp_pos[1] - y;

		double radian = atan2(p1[1], p1[0]) - atan2(p2[1], p2[0]);
		radian = abs(radian);
		if (radian > PI)
			radian = 2 * PI - radian;
		angle += radian;
	}
	if (fabs(6.28318530717958647692 - abs(angle)) < 1)
		return true;
	else
		return false;
}