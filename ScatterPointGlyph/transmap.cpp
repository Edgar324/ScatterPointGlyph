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
#include <vtkTextActor3D.h>
#include <vtkTextProperty.h>
#include <vtkPropPicker.h>
#include <QVTKWidget.h>
#include "transmap_data.h"
#include "vtkCellArray.h"
#include "tree_common.h"
#include "scatter_point_dataset.h"
#include "tour_path_generator.h"

TransMap::TransMap() {
	this->state = TransMap::NORMAL;
	this->EventCallbackCommand->SetCallback(TransMap::ProcessEvents);

	float picker_tolerence = 1e-4;
	this->level_one_node_picker = vtkPropPicker::New();
	//this->level_one_node_picker->SetTolerance(picker_tolerence);
	this->level_one_node_picker->PickFromListOn();
	
	this->level_zero_node_picker = vtkPropPicker::New();
	//this->level_zero_node_picker->SetTolerance(picker_tolerence);
	this->level_zero_node_picker->PickFromListOn();

	this->trans_picker = vtkPropPicker::New();
	//this->trans_picker->SetTolerance(picker_tolerence);
	this->trans_picker->PickFromListOn();

	this->boundary_picker = vtkPropPicker::New();
	//this->boundary_picker->SetTolerance(picker_tolerence);
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

	this->path_generator_ = new TourPathGenerator;
	this->is_map_update_needed = false;
}

TransMap::TransMap(QVTKWidget* parent) 
	: TransMap() {
	this->parent_view = parent;
}

TransMap::~TransMap() {
	
}

void TransMap::SetOriginalData(ScatterPointDataset* data) {
	this->scatter_data_ = data;
}

void TransMap::SetData(TransMapData* data) {
	dataset_ = data;

	this->path_generator_->SetData(dataset_);
	this->path_generator_->GenerateSpanningTree();
	trans_edges = this->path_generator_->edge_list;

	//this->SetEnabled(false);
	this->ResizeActors();
	this->ConstructActors();
	this->BuildRepresentation();
	this->SetEnabled(true);

	this->highlight_node_sequence.clear();
	this->UpdateHightlightActor();
	is_map_update_needed = false;
}

void TransMap::UpdateScale() {
	this->ConstructActors();
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
}

void TransMap::SetBrushSelectionOff() {
	this->state = this->state_vec[this->state_vec.size() - 1];
	this->state_vec.pop_back();
	this->selection_brush_poly->Initialize();
	this->selection_brush_actor->Modified();
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

void TransMap::SetNodeSelected(int node_id) {
	std::map< int, CNode* >::iterator iter = dataset_->cluster_node_map.begin();
	while (iter != dataset_->cluster_node_map.end() && iter->second->id != node_id) iter++;
	if (iter != dataset_->cluster_node_map.end()) {
		this->SelectNode(iter->second);
		this->UpdateHightlightActor();
	}
}

int TransMap::GetSelectedClusterIndex() {
	if (highlight_node_sequence.size() != 0) {
		std::list< CNode* >::iterator node_iter = highlight_node_sequence.begin();
		std::map< int, CNode* >::iterator iter = dataset_->cluster_node_map.begin();
		while (iter != dataset_->cluster_node_map.end() && iter->second != *node_iter) iter++;
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
		float node_center_x = node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_center_y = node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
		// construct the polydata 
		vtkPoints* points = vtkPoints::New();
		vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
		colors->SetNumberOfComponents(4);

		vtkPolyData* polydata = this->level_one_node_glyph_polys[i];
		polydata->Initialize();

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

			background_ids.push_back(points->InsertNextPoint(node_center_x + x, node_center_y + y, -0.00002));
			colors->InsertNextTuple4(128, 128, 128, 10);
		}
		polydata->InsertNextCell(VTK_POLYGON, dataset_->var_num, background_ids.data());

		// paint variance
		std::vector< int > var_point_ids;
		var_point_ids.resize(dataset_->var_num * 2);
		for (int j = 0; j < dataset_->var_num; ++j) {
			float end_arc = j * 3.14159 * 2 / dataset_->var_num;

			float temp_radius = node_radius * (node->average_values[j] - node->value_variance[j]) * 0.8;
			if (temp_radius < 0) temp_radius = 0;
			float x = temp_radius * cos(end_arc);
			float y = temp_radius * sin(end_arc);

			vtkIdType id_one = points->InsertNextPoint(node_center_x + x, node_center_y + y, -0.00001);
			colors->InsertNextTuple4(128, 128, 128, 128);

			temp_radius = node_radius * (node->average_values[j] + node->value_variance[j]) * 0.8;
			if (temp_radius > 1.0) temp_radius = 1.0;
			x = temp_radius * cos(end_arc);
			y = temp_radius * sin(end_arc);

			vtkIdType id_two = points->InsertNextPoint(node_center_x + x, node_center_y + y, -0.00001);
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
			float temp_radius = node_radius * node->average_values[j] * 0.8;
			float end_arc = j * 3.14159 * 2 / dataset_->var_num;
			float x = temp_radius * cos(end_arc);
			float y = temp_radius * sin(end_arc);

			value_pos_ids[j] = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0);
			//colors->InsertNextTuple4(dataset_->level_one_colors[3 * i], dataset_->level_one_colors[3 * i + 1], dataset_->level_one_colors[3 * i + 2], 255);
			colors->InsertNextTuple4(node->color.red(), node->color.green(), node->color.blue(), 255);
		}
		for (int j = 0; j < dataset_->var_num; ++j) {
			cell_ids[0] = value_pos_ids[j];
			cell_ids[1] = value_pos_ids[(j + 1) % dataset_->var_num];
			polydata->InsertNextCell(VTK_LINE, 2, cell_ids);
		}
		

		// paint axis
		vtkIdType center_id = points->InsertNextPoint(node_center_x, node_center_y, -0.0001);
		colors->InsertNextTuple4(200, 200, 200, 255);
		for (int j = 0; j < dataset_->var_num; ++j) {
			float end_arc = j * 3.14159 * 2 / dataset_->var_num;
			float x = node_radius * cos(end_arc) * 0.8;
			float y = node_radius * sin(end_arc) * 0.8;

			vtkIdType pre_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, -0.0001);
			colors->InsertNextTuple4(200, 200, 200, 255);
			
			cell_ids[0] = center_id;
			cell_ids[1] = pre_id;
			polydata->InsertNextCell(VTK_LINE, 2, cell_ids);
		}

		float x = node_radius * 0.9;
		float y = 0;
		float end_arc = 0;

		vtkIdType pre_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.0001);
		for (int k = 1; k <= 20; ++k) {
			end_arc = (float)k / 20 * 3.14159 * 2;
			x = node_radius * cos(end_arc) * 0.9;
			y = node_radius * sin(end_arc) * 0.9;
			vtkIdType current_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.0001);

			cell_ids[0] = pre_id;
			cell_ids[1] = current_id;
			polydata->InsertNextCell(VTK_LINE, 2, cell_ids);

			pre_id = current_id;
		}

		for (int k = 0; k < 21; ++k) colors->InsertNextTuple4(200, 200, 200, 200);

		this->level_one_node_glyph_actors[i]->Modified();
	}

	for (int i = 0; i < dataset_->level_zero_nodes.size(); ++i) {
		CNode* node = dataset_->level_zero_nodes[i];
		float node_center_x = node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_center_y = node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

		// construct the polydata
		vtkPoints* points = vtkPoints::New();
		vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
		colors->SetNumberOfComponents(3);

		vtkPolyData* polydata = this->level_zero_node_glyph_polys[i];
		polydata->Initialize();

		polydata->SetPoints(points);
		polydata->GetPointData()->SetScalars(colors);

		vtkCellArray* poly_array = vtkCellArray::New();
		polydata->SetPolys(poly_array);

		float temp_radius = node_radius * 0.2;
		float end_arc = 0;
		float x = temp_radius;
		float y = 0;

		vtkIdType center_id = points->InsertNextPoint(node_center_x, node_center_y, 0);
		vtkIdType pre_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0);

		for (int k = 1; k <= 20; ++k) {
			end_arc = (float)k / 20 * 3.14159 * 2;
			x = temp_radius * cos(end_arc);
			y = temp_radius * sin(end_arc);
			vtkIdType current_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0);

			cell_ids[0] = pre_id;
			cell_ids[1] = current_id;
			cell_ids[2] = center_id;
			polydata->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

			pre_id = current_id;
		}

		for (int k = 0; k < 22; ++k) 
			//colors->InsertNextTuple3(dataset_->level_zero_colors[3 * i], dataset_->level_zero_colors[3 * i + 1], dataset_->level_zero_colors[3 * i + 2]);
			colors->InsertNextTuple3(node->color.red(), node->color.green(), node->color.blue());

		this->level_zero_node_glyph_actors[i]->Modified();
	}

	for (int i = 0; i < trans_glyph_actors.size(); ++i){
		CNode* node_one = dataset_->level_one_nodes[trans_edges[i * 2]];
		CNode* node_two = dataset_->level_one_nodes[trans_edges[i * 2 + 1]];

		float node_one_center_x = node_one->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_one_center_y = node_one->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
		float node_two_center_x = node_two->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_two_center_y = node_two->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

		vtkPoints* points = vtkPoints::New();
		vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
		colors->SetNumberOfComponents(3);

		vtkPolyData* polydata = trans_glyph_polys[i];
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
		begin_pos[0] = node_one_center_x + node_radius * dir_vec[0];
		begin_pos[1] = node_one_center_y + node_radius * dir_vec[1];
		end_pos[0] = node_two_center_x - node_radius * dir_vec[0];
		end_pos[1] = node_two_center_y - node_radius * dir_vec[1];

		float value_dis = 0.4;

		cell_ids[0] = points->InsertNextPoint(begin_pos[0] + 0.5 * node_radius * orth_vec[0] * value_dis, begin_pos[1] + 0.5 * node_radius * orth_vec[1] * value_dis, 0);
		cell_ids[1] = points->InsertNextPoint(end_pos[0] + 0.5 * node_radius * orth_vec[0] * value_dis, end_pos[1] + 0.5 * node_radius * orth_vec[1] * value_dis, 0);
		cell_ids[2] = points->InsertNextPoint(end_pos[0], end_pos[1], 0);
		cell_ids[3] = points->InsertNextPoint(begin_pos[0], begin_pos[1], 0);

		polydata->InsertNextCell(VTK_POLYGON, 4, cell_ids);

		colors->InsertNextTuple3(163, 82, 82);
		colors->InsertNextTuple3(163, 82, 82);
		colors->InsertNextTuple3(163, 82, 82);
		colors->InsertNextTuple3(163, 82, 82);
			

		trans_glyph_actors[i]->Modified();
	}
}

void TransMap::ResizeActors() {
	if (level_one_node_glyph_actors.size() < dataset_->level_one_nodes.size()) {
		for (int i = level_one_node_glyph_actors.size(); i < dataset_->level_one_nodes.size(); ++i) {
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

			vtkActor* node_actor = vtkActor::New();
			vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
			mapper->SetInputData(polydata);
			node_actor->SetMapper(mapper);
			node_actor->GetProperty()->SetLineWidth(2.0);

			level_one_node_glyph_actors.push_back(node_actor);
			level_one_node_glyph_polys.push_back(polydata);

			if (Enabled) this->DefaultRenderer->AddActor(node_actor);

			this->level_one_node_picker->AddPickList(node_actor);
		}
	} else {
		for (int i = dataset_->level_one_nodes.size(); i < level_one_node_glyph_actors.size(); ++i) {
			if (Enabled) this->DefaultRenderer->RemoveActor(level_one_node_glyph_actors[i]);
			this->level_one_node_picker->DeletePickList(level_one_node_glyph_actors[i]);
			level_one_node_glyph_actors[i]->Delete();
		}
		level_one_node_glyph_actors.resize(dataset_->level_one_nodes.size());
		level_one_node_glyph_polys.resize(dataset_->level_one_nodes.size());
	}

	if (level_zero_node_glyph_actors.size() < dataset_->level_zero_nodes.size()) {
		for (int i = level_zero_node_glyph_actors.size(); i < dataset_->level_zero_nodes.size(); ++i) {
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

			vtkActor* node_actor = vtkActor::New();
			vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
			mapper->SetInputData(polydata);
			node_actor->SetMapper(mapper);
			node_actor->GetProperty()->SetLineWidth(2.0);

			level_zero_node_glyph_actors.push_back(node_actor);
			level_zero_node_glyph_polys.push_back(polydata);

			this->level_zero_node_picker->AddPickList(node_actor);

			if (Enabled) this->DefaultRenderer->AddActor(node_actor);
		}
	}
	else {
		for (int i = dataset_->level_zero_nodes.size(); i < level_zero_node_glyph_actors.size(); ++i) {
			if (Enabled) this->DefaultRenderer->RemoveActor(level_zero_node_glyph_actors[i]);
			this->level_zero_node_picker->DeletePickList(level_zero_node_glyph_actors[i]);
			level_zero_node_glyph_actors[i]->Delete();
		}
		level_zero_node_glyph_actors.resize(dataset_->level_zero_nodes.size());
		level_zero_node_glyph_polys.resize(dataset_->level_zero_nodes.size());
	}

	if (trans_glyph_actors.size() < path_generator_->edge_list.size() / 2) {
		for (int i = trans_glyph_actors.size(); i < path_generator_->edge_list.size() / 2; ++i) {
			vtkPoints* points = vtkPoints::New();
			vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
			colors->SetNumberOfComponents(4);

			vtkPolyData* polydata = vtkPolyData::New();
			polydata->SetPoints(points);
			polydata->GetPointData()->SetScalars(colors);

			vtkCellArray* poly_array = vtkCellArray::New();
			polydata->SetPolys(poly_array);

			vtkActor* trans_actor = vtkActor::New();
			vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
			mapper->SetInputData(polydata);
			trans_actor->SetMapper(mapper);
			trans_actor->GetProperty()->SetLineWidth(2.0);

			trans_glyph_actors.push_back(trans_actor);
			trans_glyph_polys.push_back(polydata);

			if (Enabled) this->DefaultRenderer->AddActor(trans_actor);
		}
	} else {
		for (int i = path_generator_->edge_list.size() / 2; i < trans_glyph_actors.size(); ++i) {
			if (Enabled) this->DefaultRenderer->RemoveActor(trans_glyph_actors[i]);
			trans_glyph_actors[i]->Delete();
		}
		trans_glyph_actors.resize(path_generator_->edge_list.size() / 2);
		trans_glyph_polys.resize(path_generator_->edge_list.size() / 2);
	}

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

	this->level_one_node_picker->Pick(X, Y, 1.0, this->DefaultRenderer);
	vtkProp* prop = this->level_one_node_picker->GetViewProp();
	/*vtkAssemblyPath* path = this->GetAssemblyPath(X, Y, 1.0, this->level_one_node_picker);*/
	if (prop != NULL) {
		this->state = TransMap::HIGHLIGHT_LEVEL_ONE_NODE;
		this->HighlightHandle(prop);
		return;
	}

	this->level_zero_node_picker->Pick(X, Y, 1.0, this->DefaultRenderer);
	prop = this->level_zero_node_picker->GetViewProp();
	if (prop != NULL) {
		this->state = TransMap::HIGHLIGHT_LEVEL_ZERO_NODE;
		this->HighlightHandle(prop);
		return;
	}

	this->trans_picker->Pick(X, Y, 1.0, this->DefaultRenderer);
	prop = this->trans_picker->GetViewProp();
	if (prop != NULL) {
		this->state = TransMap::HIGHLIGHT_TRANSFER_BAND;
		this->HighlightHandle(prop);
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
	if (!this->DefaultRenderer || this->state == WidgetState::SELECT_CLUSTER) {
		return;
	}

	int X = this->Interactor->GetEventPosition()[0];
	int Y = this->Interactor->GetEventPosition()[1];

	vtkAssemblyPath* path = this->GetAssemblyPath(X, Y, 0., this->level_one_node_picker);
	if (path != NULL) {
		vtkActor* actor = dynamic_cast<vtkActor*>(path->GetFirstNode()->GetViewProp());
		//if (this->current_handle == actor) return;
		if (actor == NULL) return;

		int node_index = -1;
		for (int i = 0; i < level_one_node_glyph_actors.size(); ++i)
			if (level_one_node_glyph_actors[i] == actor) {
				node_index = i;
				break;
			}
		if (node_index == -1) return;
		dataset_->level_one_nodes[node_index]->is_expanded = true;

		is_map_update_needed = true;
		return;
	}
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
	/*this->highlight_node_sequence.clear();
	for (int i = 0; i < dataset_->level_one_nodes.size(); ++i)
		if (dataset_->level_one_nodes[i]->is_highlighted)  this->highlight_node_sequence.push_back(dataset_->level_one_nodes[i]);
	for (int i = 0; i < dataset_->level_zero_nodes.size(); ++i)
		if (dataset_->level_zero_nodes[i]->is_highlighted) this->highlight_node_sequence.push_back(dataset_->level_zero_nodes[i]);*/
	if (seqence_text_actors.size() < highlight_node_sequence.size()) {
		for (int i = 0; i < highlight_node_sequence.size() - seqence_text_actors.size(); ++i) {
			vtkTextActor3D* actor = vtkTextActor3D::New();
			actor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
			actor->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
			actor->GetTextProperty()->SetBold(true);
			actor->GetTextProperty()->SetFontFamilyToArial();
			//actor->GetTextProperty()->SetFontSize(20);
			//actor->SetPosition(0, 0, 0.01);
			//actor->SetInput("Test");
			//actor->Modified();
			this->DefaultRenderer->AddActor(actor);
			seqence_text_actors.push_back(actor);
		}
	} else {
		for (int i = seqence_text_actors.size(); i >= highlight_node_sequence.size() + 1; --i) {
			this->DefaultRenderer->RemoveActor(seqence_text_actors[i - 1]);
			seqence_text_actors[i - 1]->Delete();
		}
		seqence_text_actors.resize(highlight_node_sequence.size());
	}

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
	int hightlight_node_index = 0;
	while (iter != highlight_node_sequence.end()) {
		int node_index = -1;
		for (int i = 0; i < dataset_->level_one_nodes.size(); ++i)
			if (dataset_->level_one_nodes[i] == *iter) {
				node_index = i;
				break;
			}
		if (node_index != -1) {
			CNode* temp_node = dataset_->level_one_nodes[node_index];
			float node_center_x = temp_node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
			float node_center_y = temp_node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

			vtkIdType pre_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.0001);
			for (int k = 1; k <= 20; ++k) {
				end_arc = (float)k / 20 * 3.14159 * 2;
				x = node_radius * cos(end_arc);
				y = node_radius * sin(end_arc);
				vtkIdType current_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.0001);

				cell_ids[0] = pre_id;
				cell_ids[1] = current_id;
				this->highlight_poly->InsertNextCell(VTK_LINE, 2, cell_ids);

				pre_id = current_id;
			}

			for (int k = 0; k < 21; ++k) color_array->InsertNextTuple3(255, 0, 0);
			iter++;

			char buffer[20];
			itoa(hightlight_node_index, buffer, 10);
			seqence_text_actors[hightlight_node_index]->SetInput(buffer);
			seqence_text_actors[hightlight_node_index]->SetPosition(node_center_x - node_radius, node_center_y + node_radius * 0.5, 0.001);
			seqence_text_actors[hightlight_node_index]->GetTextProperty()->SetFontSize(20);
			float scale = node_radius * 0.5 / 20;
			seqence_text_actors[hightlight_node_index]->SetScale(scale, scale, scale);
			seqence_text_actors[hightlight_node_index]->Modified();
			hightlight_node_index++;

			continue;
		}

		for (int i = 0; i < dataset_->level_zero_nodes.size(); ++i)
			if (dataset_->level_zero_nodes[i] == *iter) {
				node_index = i;
				break;
			}
		if (node_index != -1) {
			CNode* temp_node = dataset_->level_zero_nodes[node_index];
			float node_center_x = temp_node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
			float node_center_y = temp_node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

			vtkIdType pre_id = points->InsertNextPoint(node_center_x + x * 0.4, node_center_y + y * 0.4, 0.0001);
			for (int k = 1; k <= 20; ++k) {
				end_arc = (float)k / 20 * 3.14159 * 2;
				x = node_radius * cos(end_arc);
				y = node_radius * sin(end_arc);
				vtkIdType current_id = points->InsertNextPoint(node_center_x + x * 0.4, node_center_y + y * 0.4, 0.0001);

				cell_ids[0] = pre_id;
				cell_ids[1] = current_id;
				this->highlight_poly->InsertNextCell(VTK_LINE, 2, cell_ids);

				pre_id = current_id;
			}

			for (int k = 0; k < 21; ++k) color_array->InsertNextTuple3(255, 0, 0);
			iter++;

			char buffer[20];
			itoa(hightlight_node_index, buffer, 10);
			seqence_text_actors[hightlight_node_index]->SetInput(buffer);
			seqence_text_actors[hightlight_node_index]->SetPosition(node_center_x - node_radius, node_center_y + node_radius * 0.5, 0.001);
			seqence_text_actors[hightlight_node_index]->GetTextProperty()->SetFontSize(20);
			float scale = node_radius * 0.5 / 20;
			seqence_text_actors[hightlight_node_index]->SetScale(scale, scale, scale);
			seqence_text_actors[hightlight_node_index]->Modified();
			hightlight_node_index++;

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
			std::list< CNode* >::iterator iter = highlight_node_sequence.begin();
			while (iter != highlight_node_sequence.end()) {
				(*iter)->is_highlighted = false;
				iter++;
			}
			highlight_node_sequence.clear();
			highlight_node_sequence.push_back(node);
			node->is_highlighted = true;
		} else if (highlight_node_sequence.size() == 1) {
			if (*highlight_node_sequence.begin() == node) {
				std::list< CNode* >::iterator iter = highlight_node_sequence.begin();
				while (iter != highlight_node_sequence.end()) {
					(*iter)->is_highlighted = false;
					iter++;
				}
				highlight_node_sequence.clear();
			}
			else {
				std::list< CNode* >::iterator iter = highlight_node_sequence.begin();
				while (iter != highlight_node_sequence.end()) {
					(*iter)->is_highlighted = false;
					iter++;
				}
				highlight_node_sequence.clear();
				highlight_node_sequence.push_back(node);
				node->is_highlighted = true;
			}
		} else {
			highlight_node_sequence.push_back(node);
			node->is_highlighted = true;
		}
	} else {
		std::list< CNode* >::iterator iter = highlight_node_sequence.begin();
		while (iter != highlight_node_sequence.end() && *iter != node) iter++;
		if (iter != highlight_node_sequence.end()) {
			(*iter)->is_highlighted = false;
			highlight_node_sequence.erase(iter);
		} else {
			highlight_node_sequence.push_back(node);
			node->is_highlighted = true;
		}
	}
}

void TransMap::UpdateBrushSelectionResult() {
	std::list< CNode* >::iterator iter = highlight_node_sequence.begin();
	while (iter != highlight_node_sequence.end()) {
		(*iter)->is_highlighted = false;
		iter++;
	}
	highlight_node_sequence.clear();

	for (int i = 0; i < dataset_->level_one_nodes.size(); ++i) {
		CNode* node = dataset_->level_one_nodes[i];
		float node_center_x = node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_center_y = node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
		if (IsInsideSelection(node_center_x, node_center_y)) {
			highlight_node_sequence.push_back(dataset_->level_one_nodes[i]);
			dataset_->level_one_nodes[i]->is_highlighted = true;
		}
	}

	for (int i = 0; i < dataset_->level_zero_nodes.size(); ++i) {
		CNode* node = dataset_->level_zero_nodes[i];
		float node_center_x = node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_center_y = node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
		if (IsInsideSelection(node_center_x, node_center_y)) {
			highlight_node_sequence.push_back(dataset_->level_zero_nodes[i]);
			dataset_->level_zero_nodes[i]->is_highlighted = true;
		}
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