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
#include <vtkTooltipItem.h>
#include <QVTKWidget.h>
#include "transmap_data.h"
#include "vtkCellArray.h"
#include "tree_common.h"
#include "scatter_point_dataset.h"
#include "tour_path_generator.h"
#include "scatter_point_view.h"

//#define RADAR_GLYPH

TransMap::TransMap() {
	this->state = TransMap::SELECT_SINGLE_CLUSTER;
	this->EventCallbackCommand->SetCallback(TransMap::ProcessEvents);

	float picker_tolerence = 1e-4;
	this->node_picker = vtkPropPicker::New();
	//this->level_one_node_picker->SetTolerance(picker_tolerence);
	this->node_picker->PickFromListOn();

	double bounds[6];
	bounds[0] = -0.5;
	bounds[1] = 0.5;
	bounds[2] = -0.5;
	bounds[3] = 0.5;
	bounds[4] = -0.5;
	bounds[5] = 0.5;
	this->PlaceFactor = 1.0;
	this->PlaceWidget(bounds);

	this->node_radius_ = 10;
	this->current_node_ = NULL;
	this->scatter_data_ = NULL;
	this->dataset_ = NULL;
    this->glyph_indicator_renderer_ = NULL;

	this->highlight_actor = vtkActor::New();
	this->highlight_poly = vtkPolyData::New();
	this->hightlight_mapper = vtkPolyDataMapper::New();
	this->hightlight_mapper->SetInputData(this->highlight_poly);
	this->highlight_actor->SetMapper(this->hightlight_mapper);
	this->highlight_actor->GetProperty()->SetLineWidth(3.0);
	
	this->trans_glyph_actors = vtkActor::New();
	this->trans_glyph_polys = vtkPolyData::New();
	this->trans_glyph_mapper = vtkPolyDataMapper::New();
	this->trans_glyph_mapper->SetInputData(this->trans_glyph_polys);
	this->trans_glyph_actors->SetMapper(this->trans_glyph_mapper);
	this->trans_glyph_actors->GetProperty()->SetColor(1.0, 0.55, 0.24);
    //this->trans_glyph_actors->GetProperty()->SetColor(0.6, 0.6, 0.6);

    this->linked_glyph_actors = vtkActor::New();
	this->linked_glyph_polys = vtkPolyData::New();
	this->linked_glyph_mapper = vtkPolyDataMapper::New();
	this->linked_glyph_mapper->SetInputData(this->linked_glyph_polys);
	this->linked_glyph_actors->SetMapper(this->linked_glyph_mapper);
	this->linked_glyph_actors->GetProperty()->SetColor(0.15, 0.65, 0.3);

	this->selection_brush_poly = vtkPolyData::New();
	this->selection_brush_mapper = vtkPolyDataMapper::New();
	this->selection_brush_actor = vtkActor::New();
	this->selection_brush_mapper->SetInputData(this->selection_brush_poly);
	this->selection_brush_actor->SetMapper(this->selection_brush_mapper);
	this->selection_brush_actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
	this->selection_brush_actor->GetProperty()->SetLineWidth(3.0);

    this->density_actor_ = vtkActor::New();
	this->density_poly_ = vtkPolyData::New();
	this->density_mapper_ = vtkPolyDataMapper::New();
	this->density_mapper_->SetInputData(this->density_poly_);
	this->density_actor_->SetMapper(this->density_mapper_);
    density_actor_->SetVisibility(false);

    this->indicator_actor_ = vtkActor::New();
	this->indicator_poly_ = vtkPolyData::New();
	this->indicator_mapper_ = vtkPolyDataMapper::New();
	this->indicator_mapper_->SetInputData(this->indicator_poly_);
	this->indicator_actor_->SetMapper(this->indicator_mapper_);
    this->indicator_actor_->GetProperty()->SetColor(0.7, 0.7, 0.7);

    this->var_icon_actor_ = vtkActor::New();
	this->var_icon_poly_ = vtkPolyData::New();
	this->var_icon_mapper_ = vtkPolyDataMapper::New();
	this->var_icon_mapper_->SetInputData(this->var_icon_poly_);
	this->var_icon_actor_->SetMapper(this->var_icon_mapper_);

    indicator_text_ = vtkTextActor3D::New();
    indicator_text_->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
	indicator_text_->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
	indicator_text_->GetTextProperty()->SetBold(true);
	indicator_text_->GetTextProperty()->SetFontFamilyToArial();

	//this->tool_tip_item_ = vtkTooltipItem::New();

	this->path_generator_ = new TourPathGenerator;
}

TransMap::TransMap(ScatterPointView* parent)
	: TransMap() {
	this->parent_view = parent;
}

TransMap::~TransMap() {
	
}

void TransMap::SetData(ScatterPointDataset* ori_data, TransMapData* data) {
	this->scatter_data_ = ori_data;
	this->dataset_ = data;

	axis_order_.resize(this->scatter_data_->var_num);
	for (int i = 0; i < this->scatter_data_->var_num; ++i) axis_order_[i] = i;

	this->path_generator_->SetData(dataset_);

	this->highlight_node_sequence.clear();
	this->current_node_ = NULL;
	this->trans_edges.clear();

	if (is_mst_fixed_) this->path_generator_->GenerateSpanningTree();
	if (is_var_trend_fixed_ && var_trend_index_ != -1) {
		this->path_generator_->GenerateVarTrend(var_trend_index_);

		this->highlight_node_sequence.clear();
		for (int i = 0; i < this->path_generator_->trans_edge_list.size() / 2; ++i) {
			this->highlight_node_sequence.push_back(dataset_->level_one_nodes[this->path_generator_->trans_edge_list[i * 2]]);
		}
		if (this->path_generator_->trans_edge_list.size() != 0) {
			this->highlight_node_sequence.push_back(dataset_->level_one_nodes[this->path_generator_->trans_edge_list[this->path_generator_->trans_edge_list.size() - 1]]);
		}
	}

	this->BuildRepresentation();
	this->SetEnabled(true);
}

void TransMap::SetNodeRadius(float r) {
	this->node_radius_ = r;
}

void TransMap::ShowMinimumSpanningTree(bool enabled) {
	this->linked_edges.clear();
	this->is_mst_fixed_ = enabled;
	if (enabled) this->path_generator_->GenerateSpanningTree();
	this->UpdateTransEdgeActor();
}

void TransMap::ShowVarTrend(int var_index) {
	this->trans_edges.clear();
	if (var_index < 0)
		this->is_var_trend_fixed_ = false;
	else
		this->is_var_trend_fixed_ = true;

	var_trend_index_ = var_index;
	if (var_index >= 0) this->path_generator_->GenerateVarTrend(var_index);

	this->highlight_node_sequence.clear();
	if (this->is_var_trend_fixed_) {
		for (int i = 0; i < this->path_generator_->trans_edge_list.size() / 2; ++i) {
			this->highlight_node_sequence.push_back(dataset_->level_one_nodes[this->path_generator_->trans_edge_list[i * 2]]);
		}
		if (this->path_generator_->trans_edge_list.size() != 0) {
			this->highlight_node_sequence.push_back(dataset_->level_one_nodes[this->path_generator_->trans_edge_list[this->path_generator_->trans_edge_list.size() - 1]]);
		}
	}

	this->UpdateHightlightActor();
	this->UpdateTransEdgeActor();
}

void TransMap::HighlightVar(int var_index) {
    bool is_exist = false;
    for (int i = 0; i < focus_var_index_.size(); ++i) {
        if (focus_var_index_[i] == var_index) is_exist = true;
    }
    if (!is_exist) focus_var_index_.push_back(var_index);
	if (this->focus_var_index_.size() != 0) {
		this->is_focus_var_fixed_ = true;
	} else {
		this->is_focus_var_fixed_ = false;
	}

	this->UpdateNodeActors();
	this->parent_view->update();
}

void TransMap::SetInteractionState(WidgetState s){
	this->state = s;

	while (this->highlight_node_sequence.size() > 0) {
		this->highlight_node_sequence.front()->is_highlighted = false;
		this->highlight_node_sequence.pop_front();
	}
	this->current_node_ = NULL;
	this->UpdateHightlightActor();

	this->trans_edges.clear();
	this->UpdateTransEdgeActor();

	this->parent_view->update();
}

void TransMap::SetAxisOrder(std::vector< int >& order) {
	axis_order_ = order;
	this->UpdateNodeActors();
}

void TransMap::SetIndicatorRenderer(vtkRenderer* renderer) {
    glyph_indicator_renderer_ = renderer;
}

void TransMap::OnNodeSelected(int node_id) {
	std::map< int, CNode* >::iterator iter = dataset_->cluster_node_map.begin();
    while (iter != dataset_->cluster_node_map.end() && iter->second->id() != node_id) iter++; 
	if (iter != dataset_->cluster_node_map.end()) OnNodeSelected(iter->second);
}

void TransMap::GetFocusVarIndex(std::vector< int >& focus_index) {
    focus_index = this->focus_var_index_;
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

void TransMap::GetSelectedClusterIds(std::vector< int >& ids) {
	ids.clear();

	std::list< CNode* >::iterator node_iter = highlight_node_sequence.begin();
	while (node_iter != highlight_node_sequence.end()) {
		ids.push_back((*node_iter)->id());
		node_iter++;
	}
}

void TransMap::GetSelectedClusterNodes(std::vector< CNode* >& nodes) {
	nodes.clear();

	std::list< CNode* >::iterator node_iter = highlight_node_sequence.begin();
	while (node_iter != highlight_node_sequence.end()) {
		nodes.push_back(*node_iter);
		node_iter++;
	}
}

void TransMap::UpdateNodeActors() {
	// resize node actors
	if (node_glyph_actors.size() < dataset_->cluster_nodes.size()) {
		for (int i = node_glyph_actors.size(); i < dataset_->cluster_nodes.size(); ++i) {
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
			node_actor->GetProperty()->SetLineWidth(3.0);
			//node_actor->GetProperty()->SetLineStipplePattern(0xF00F);

			node_glyph_actors.push_back(node_actor);
			node_glyph_polys.push_back(polydata);

			if (Enabled) this->DefaultRenderer->AddActor(node_actor);

			this->node_picker->AddPickList(node_actor);

			vtkTextActor3D* actor = vtkTextActor3D::New();
			actor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
			actor->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
			actor->GetTextProperty()->SetBold(true);
			actor->GetTextProperty()->SetFontFamilyToArial();
			this->DefaultRenderer->AddActor(actor);
		}
	}
	else {
		for (int i = dataset_->cluster_nodes.size(); i < node_glyph_actors.size(); ++i) {
			if (Enabled) {
				this->DefaultRenderer->RemoveActor(node_glyph_actors[i]);
			}
			this->node_picker->DeletePickList(node_glyph_actors[i]);
			node_glyph_actors[i]->Delete();
		}
		node_glyph_actors.resize(dataset_->cluster_nodes.size());
		node_glyph_polys.resize(dataset_->cluster_nodes.size());
	}

    int max_node_point_num = -1;
    for (int i = 0; i < dataset_->level_one_nodes.size(); ++i)
        if (dataset_->level_one_nodes[i]->point_count > max_node_point_num)
            max_node_point_num = dataset_->level_one_nodes[i]->point_count;

    // update node actors
	vtkIdType cell_ids[5];

    // update indicator actor
    vtkPoints* indicator_points = vtkPoints::New();
	indicator_poly_->Initialize();
	indicator_poly_->SetPoints(indicator_points);
	vtkCellArray* indicator_poly_array = vtkCellArray::New();
	indicator_poly_->SetPolys(indicator_poly_array);

    int seg_per_circle = 50;
    int seg_num = seg_per_circle - 1;
    float gray_value = 128;
    std::vector< vtkIdType > inner_ids, outer_ids;
    for (int j = 0; j <= seg_num; ++j) {
        float end_arc = -1 * j * 3.14159 * 2 / seg_per_circle + 3.14159 * 0.5;
        float x = 1.0 * cos(end_arc) * 1.05;
		float y = 1.0 * sin(end_arc) * 1.05;

        vtkIdType id_one = indicator_points->InsertNextPoint(x, y, 0.000);
        inner_ids.push_back(id_one);

        x = 1.0 * cos(end_arc) * 1.18;
		y = 1.0 * sin(end_arc) * 1.18;

        vtkIdType id_two = indicator_points->InsertNextPoint(x, y, 0.000);
        outer_ids.push_back(id_two);
    }

    for (int j = 0; j < outer_ids.size() - 1; ++j) {
        cell_ids[0] = outer_ids[j];
        cell_ids[1] = outer_ids[j + 1];
        cell_ids[2] = inner_ids[j];
        indicator_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

        cell_ids[0] = outer_ids[j + 1];
        cell_ids[1] = inner_ids[j + 1];
        cell_ids[2] = inner_ids[j];
        indicator_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
    }
    indicator_actor_->Modified();

    char buffer[20];
	itoa(max_node_point_num, buffer, 10);
    std::string str = std::string("Max Point Num: ") + std::string(buffer);
	indicator_text_->SetInput(str.c_str());
	indicator_text_->SetPosition(-1, -1.5, 0);
	indicator_text_->GetTextProperty()->SetFontSize(20);
	float scale = 0.2 / 20;
	indicator_text_->SetScale(scale, scale, scale);
	indicator_text_->Modified();
    
	for (int i = 0; i < dataset_->cluster_nodes.size(); ++i) {
		CNode* node = dataset_->cluster_nodes[i];
		float node_center_x = node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_center_y = node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

		if (node->point_count >= dataset_->min_point_num) {
			// construct the polydata 
			vtkPoints* points = vtkPoints::New();
			vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
			colors->SetNumberOfComponents(4);

			vtkPolyData* polydata = this->node_glyph_polys[i];
			polydata->Initialize();

			polydata->SetPoints(points);
			polydata->GetPointData()->SetScalars(colors);

			vtkCellArray* poly_array = vtkCellArray::New();
			polydata->SetPolys(poly_array);

			vtkCellArray* line_array = vtkCellArray::New();
			polydata->SetLines(line_array);

			vtkCellArray* strip_array = vtkCellArray::New();
			polydata->SetStrips(strip_array);

            int seg_per_circle = 50;
			// insert background
			std::vector< vtkIdType > background_ids;
			for (int j = 0; j <= seg_per_circle; ++j) {
				float end_arc = j * 3.14159 * 2 / seg_per_circle;
				float x = node_radius_ * cos(end_arc);
				float y = node_radius_ * sin(end_arc);

				background_ids.push_back(points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.000));
				
                colors->InsertNextTuple4(255, 255, 255, 255);
                /*if (is_densitymap_shown_)
                    colors->InsertNextTuple4(255, 255, 255, 255);
                else
                    colors->InsertNextTuple4(255, 255, 255, 10);*/
			}
			polydata->InsertNextCell(VTK_POLYGON, seg_per_circle + 1, background_ids.data());

            // insert background contour
            std::vector< vtkIdType > background_contour_ids;
			for (int j = 0; j <= seg_per_circle; ++j) {
				float end_arc = j * 3.14159 * 2 / seg_per_circle;
				float x = node_radius_ * cos(end_arc);
				float y = node_radius_ * sin(end_arc);

				background_contour_ids.push_back(points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.000));
				
                colors->InsertNextTuple4(200, 200, 200, 255);
			}
            for (int j = 0; j < background_contour_ids.size() - 1; ++j) {
                cell_ids[0] = background_contour_ids[j];
                cell_ids[1] = background_contour_ids[j + 1];
                polydata->InsertNextCell(VTK_LINE, 2, cell_ids);
            }

			std::vector< vtkIdType > center_cirlce_ids;
			for (int j = 0; j <= seg_per_circle; ++j) {
				float end_arc = j * 3.14159 * 2 / seg_per_circle;
				float x = node_radius_ * 0.05 * cos(end_arc);
				float y = node_radius_ * 0.05 * sin(end_arc);

				center_cirlce_ids.push_back(points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.003));
				colors->InsertNextTuple4(0, 0, 0, 255);
			}
			polydata->InsertNextCell(VTK_POLYGON, seg_per_circle + 1, center_cirlce_ids.data());

            // paint saliency circle
			/*std::vector< vtkIdType > radius_circle_ids;
            int segment_per_circle_temp = 20;
			for (int j = 0; j <= segment_per_circle_temp; ++j) {
				float end_arc = j * 3.14159 * 2 / segment_per_circle_temp;
				float x = node_radius_ * cos(end_arc) * 1.03;
				float y = node_radius_ * sin(end_arc) * 1.03;

				radius_circle_ids.push_back(points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.0001));
                float gray_value = 250 * (1.0 - exp(-3 * (1.0 - node->saliency)));
				colors->InsertNextTuple4(gray_value, gray_value, gray_value, 255);
			}

			vtkIdType circle_ids[2];
			for (int j = 0; j < segment_per_circle_temp; ++j) {
				circle_ids[0] = radius_circle_ids[j];
				circle_ids[1] = radius_circle_ids[j + 1];
				polydata->InsertNextCell(VTK_LINE, 2, circle_ids);
			}*/

            // paint encoding for the point number
            float point_rate = (float)node->point_count / max_node_point_num;
            int seg_num = (seg_per_circle - 1) * point_rate;
            float gray_value = 250 * (1.0 - exp(-3 * (1.0 - node->saliency)));
            std::vector< vtkIdType > inner_ids, outer_ids;
            for (int j = 0; j <= seg_num; ++j) {
                float end_arc = -1 * j * 3.14159 * 2 / seg_per_circle + 3.14159 * 0.5;
                float x = node_radius_ * cos(end_arc) * 1.05;
				float y = node_radius_ * sin(end_arc) * 1.05;

                vtkIdType id_one = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
				colors->InsertNextTuple4(gray_value, gray_value, gray_value, 255);
                inner_ids.push_back(id_one);

                x = node_radius_ * cos(end_arc) * 1.18;
				y = node_radius_ * sin(end_arc) * 1.18;

                vtkIdType id_two = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
				colors->InsertNextTuple4(gray_value, gray_value, gray_value, 255);
                outer_ids.push_back(id_two);
            }
            for (int j = 0; j < outer_ids.size() - 1; ++j) {
                cell_ids[0] = outer_ids[j];
                cell_ids[1] = outer_ids[j + 1];
                cell_ids[2] = inner_ids[j];
                polydata->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

                cell_ids[0] = outer_ids[j + 1];
                cell_ids[1] = inner_ids[j + 1];
                cell_ids[2] = inner_ids[j];
                polydata->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
            }

#ifdef RADAR_GLYPH
			// paint variance
			std::vector< int > var_point_ids;
			var_point_ids.resize(dataset_->var_num * 2);
			for (int j = 0; j < dataset_->var_num; ++j) {
                int var_index = axis_order_[j];
				float end_arc = j * 3.14159 * 2 / dataset_->var_num;

				float temp_radius = node_radius_ * (node->average_values[var_index] - node->variable_variances[var_index]) * 0.9 + node_radius_ * 0.1;
				if (temp_radius < 0) temp_radius = 0;
				float x = temp_radius * cos(end_arc);
				float y = temp_radius * sin(end_arc);

				vtkIdType id_one = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
				colors->InsertNextTuple4(128, 128, 128, 128);

				temp_radius = node_radius_ * (node->average_values[var_index] + node->variable_variances[var_index]) * 0.9 + node_radius_ * 0.1;
				if (temp_radius > node_radius_) temp_radius = node_radius_;
				x = temp_radius * cos(end_arc);
				y = temp_radius * sin(end_arc);

				vtkIdType id_two = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
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
                int var_index = axis_order_[j];
				float temp_radius = node_radius_ * node->average_values[var_index] * 0.9 + node_radius_ * 0.1;
				float end_arc = j * 3.14159 * 2 / dataset_->var_num;
				float x = temp_radius * cos(end_arc);
				float y = temp_radius * sin(end_arc);

				value_pos_ids[j] = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
				//colors->InsertNextTuple4(dataset_->level_one_colors[3 * i], dataset_->level_one_colors[3 * i + 1], dataset_->level_one_colors[3 * i + 2], 255);
				colors->InsertNextTuple4(node->color.red(), node->color.green(), node->color.blue(), 255);
			}
			for (int j = 0; j < dataset_->var_num; ++j) {
				cell_ids[0] = value_pos_ids[j];
				cell_ids[1] = value_pos_ids[(j + 1) % dataset_->var_num];
				polydata->InsertNextCell(VTK_LINE, 2, cell_ids);
			}


			// paint axis
			vtkIdType center_id = points->InsertNextPoint(node_center_x, node_center_y, 0.0015);
			colors->InsertNextTuple4(200, 200, 200, 255);
			for (int j = 0; j < dataset_->var_num; ++j) {
				float end_arc = j * 3.14159 * 2 / dataset_->var_num;
				float x = node_radius_ * cos(end_arc);
				float y = node_radius_ * sin(end_arc);

				vtkIdType pre_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.0015);
				colors->InsertNextTuple4(200, 200, 200, 255);

				cell_ids[0] = center_id;
				cell_ids[1] = pre_id;
				polydata->InsertNextCell(VTK_LINE, 2, cell_ids);
			}
#else
			// paint variance
			std::vector< vtkIdType > var_point_ids1;
			std::vector< vtkIdType > var_point_ids2;
			std::vector< vtkIdType > gray_region_ids;
			int seg_per_pie = 30;
			var_point_ids1.resize(2 * seg_per_pie);
			var_point_ids2.resize(2 * seg_per_pie);
			gray_region_ids.resize(seg_per_pie + 1);

			vtkIdType line_ids[2];
			float alpha = 1.0;

			for (int j = 0; j < dataset_->var_num; ++j) {
				float begin_arc = j * 3.14159 * 2 / dataset_->var_num;
				float end_arc = (j + 1) * 3.14159 * 2 / dataset_->var_num;
				float step_arc = (end_arc - begin_arc) / (seg_per_pie - 1);

				int var_index = axis_order_[j];
                QColor pie_color = dataset_->dataset->var_colors[var_index];

                bool is_var_highlighted = false;
                for (int k = 0; k < focus_var_index_.size(); ++k)
                    if (focus_var_index_[k] == var_index) {
                        is_var_highlighted = true;
                        break;
                    }
				if (is_var_highlighted) alpha = 1.0;
				else {
                    if (focus_var_index_.size() != 0) {
                        if (is_hovering_) alpha = 0.5;
                        else alpha = 0.0;
                    }
                    else {
                        alpha = 1.0;
                    }
				}

				float temp_arc = begin_arc;
				for (int k = 0; k < seg_per_pie; ++k) {
					float cos_value = cos(temp_arc);
					float sin_value = sin(temp_arc);

					// insert the outer polygon
					float temp_radius = node_radius_ * (node->average_values[var_index] + node->variable_variances[var_index]) * 0.9 + node_radius_ * 0.1;
					if (temp_radius > node_radius_) temp_radius = node_radius_;
					float x = temp_radius * cos_value;
					float y = temp_radius * sin_value;
					vtkIdType id_one = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
					colors->InsertNextTuple4(pie_color.red(), pie_color.green(), pie_color.blue(), 128 * alpha);
					var_point_ids1[2 * k + 1] = id_one;

					// insert the inner polygon
					temp_radius = node_radius_ * (node->average_values[var_index] - node->variable_variances[var_index]) * 0.9 + node_radius_ * 0.1;
					if (temp_radius < 0) temp_radius = 0;
					x = temp_radius * cos_value;
					y = temp_radius * sin_value;

					vtkIdType id_two = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
					colors->InsertNextTuple4(pie_color.red(), pie_color.green(), pie_color.blue(), 128 * alpha);
					var_point_ids2[2 * k + 1] = id_two;

					// insert the gray polygon
					vtkIdType id_four = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
					colors->InsertNextTuple4(230, 230, 230, 255 * alpha);
					gray_region_ids[k] = id_four;

					// insert the average value
					temp_radius = node_radius_ * node->average_values[var_index] * 0.9 + node_radius_ * 0.1;
					x = temp_radius * cos_value;
					y = temp_radius * sin_value;
					vtkIdType id_three = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
					colors->InsertNextTuple4(pie_color.red(), pie_color.green(), pie_color.blue(), 255 * alpha);
					var_point_ids1[2 * k] = id_three;
					var_point_ids2[2 * k] = id_three;

					temp_arc += step_arc;
				}
				
				vtkIdType center_id = points->InsertNextPoint(node_center_x, node_center_y, 0.002);
				colors->InsertNextTuple4(230, 230, 230, 255 * alpha);
				gray_region_ids[seg_per_pie] = center_id;

				polydata->InsertNextCell(VTK_POLYGON, (int)gray_region_ids.size(), gray_region_ids.data());
				polydata->InsertNextCell(VTK_TRIANGLE_STRIP, (int)var_point_ids1.size(), var_point_ids1.data());
				polydata->InsertNextCell(VTK_TRIANGLE_STRIP, (int)var_point_ids2.size(), var_point_ids2.data());

				/*for (int k = 0; k < seg_per_pie - 1; ++k) {
					line_ids[0] = var_point_ids1[2 * k];
					line_ids[1] = var_point_ids1[2 * (k + 1)];
					polydata->InsertNextCell(VTK_LINE, 2, line_ids);
				}*/

				if (var_index == current_highlight_var_index_ && current_node_ != NULL) {
					std::vector< int > comp_ids;
					comp_ids.resize(seg_per_pie);

					float temp_arc = begin_arc;
					for (int k = 0; k < seg_per_pie; ++k) {
						float cos_value = cos(temp_arc);
						float sin_value = sin(temp_arc);

						// insert the average value
						float temp_radius = node_radius_ * current_node_->average_values[var_index] * 0.8 + node_radius_ * 0.1;
						float x = temp_radius * cos_value;
						float y = temp_radius * sin_value;
						comp_ids[k] = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.0025);
						colors->InsertNextTuple4(255, 0, 0, 255);

						temp_arc += step_arc;
					}

					for (int k = 0; k < seg_per_pie - 1; ++k) {
						line_ids[0] = comp_ids[k];
						line_ids[1] = comp_ids[k + 1];
						polydata->InsertNextCell(VTK_LINE, 2, line_ids);
					}
				}
			}
#endif // RADAR_GLYPH
		} else {
			vtkPoints* points = vtkPoints::New();
			vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
			colors->SetNumberOfComponents(3);

			vtkPolyData* polydata = this->node_glyph_polys[i];
			polydata->Initialize();

			polydata->SetPoints(points);
			polydata->GetPointData()->SetScalars(colors);

			vtkCellArray* poly_array = vtkCellArray::New();
			polydata->SetPolys(poly_array);

			float temp_radius = node_radius_ * 0.2;
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
				colors->InsertNextTuple3(node->color.red(), node->color.green(), node->color.blue());
		}

		this->node_glyph_actors[i]->Modified();
	}
}

void TransMap::UpdateTransEdgeActor() {
    this->linked_edges.clear();

    vtkIdType cell_ids[5];

	if (is_mst_fixed_) {
		for (int i = 0; i < this->path_generator_->linked_edge_list.size(); ++i) {
			int node_index = GetClusterNodeIndex(dataset_->level_one_nodes[this->path_generator_->linked_edge_list[i]]);
			this->linked_edges.push_back(node_index);
		}
    }

    this->linked_glyph_polys->Initialize();

	vtkPoints* linked_points = vtkPoints::New();
	linked_glyph_polys->SetPoints(linked_points);

	vtkCellArray* linked_poly_array = vtkCellArray::New();
	linked_glyph_polys->SetPolys(linked_poly_array);

	for (int i = 0; i < linked_edges.size() / 2; ++i){
		CNode* node_one = dataset_->cluster_nodes[linked_edges[i * 2]];
		CNode* node_two = dataset_->cluster_nodes[linked_edges[i * 2 + 1]];

		float node_one_center_x = node_one->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_one_center_y = node_one->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
		float node_two_center_x = node_two->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_two_center_y = node_two->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

		float dir_vec[2], orth_vec[2];
		dir_vec[0] = node_two->center_pos[0] - node_one->center_pos[0];
		dir_vec[1] = node_two->center_pos[1] - node_one->center_pos[1];
		float dir_length = sqrt(pow(dir_vec[0], 2) + pow(dir_vec[1], 2));
		dir_vec[0] /= dir_length;
		dir_vec[1] /= dir_length;
		orth_vec[0] = -1 * dir_vec[1];
		orth_vec[1] = dir_vec[0];

		float begin_pos[2], end_pos[2];
		begin_pos[0] = node_one_center_x + node_radius_ * dir_vec[0] * 1.2;
		begin_pos[1] = node_one_center_y + node_radius_ * dir_vec[1] * 1.2;
		end_pos[0] = node_two_center_x - node_radius_ * dir_vec[0] * 1.2;
		end_pos[1] = node_two_center_y - node_radius_ * dir_vec[1] * 1.2;

		float value_dis = 0.0;
        for (int j = 0; j < scatter_data_->var_num; ++j)
            value_dis += abs(node_one->average_values[j] - node_two->average_values[j]) * scatter_data_->var_weights[j];

		cell_ids[0] = linked_points->InsertNextPoint(begin_pos[0] + node_radius_ * orth_vec[0] * value_dis, begin_pos[1] + node_radius_ * orth_vec[1] * value_dis, 0);
		cell_ids[1] = linked_points->InsertNextPoint(end_pos[0] + node_radius_ * orth_vec[0] * value_dis, end_pos[1] + node_radius_ * orth_vec[1] * value_dis, 0);
		cell_ids[2] = linked_points->InsertNextPoint(end_pos[0] - node_radius_ * orth_vec[0] * value_dis, end_pos[1] - node_radius_ * orth_vec[1] * value_dis, 0);
		cell_ids[3] = linked_points->InsertNextPoint(begin_pos[0] - node_radius_ * orth_vec[0] * value_dis, begin_pos[1] - node_radius_ * orth_vec[1] * value_dis, 0);

		linked_glyph_polys->InsertNextCell(VTK_POLYGON, 4, cell_ids);
	}

	this->linked_glyph_actors->Modified();

    if (is_var_trend_fixed_) {
        this->trans_edges.clear();
		for (int i = 0; i < this->path_generator_->trans_edge_list.size(); ++i) {
			int node_index = GetClusterNodeIndex(dataset_->level_one_nodes[this->path_generator_->trans_edge_list[i]]);
			this->trans_edges.push_back(node_index);
		}
    }

	this->trans_glyph_polys->Initialize();

	vtkPoints* trans_points = vtkPoints::New();
	trans_glyph_polys->SetPoints(trans_points);

	vtkCellArray* trans_poly_array = vtkCellArray::New();
	trans_glyph_polys->SetPolys(trans_poly_array);

	for (int i = 0; i < trans_edges.size() / 2; ++i){
		CNode* node_one = dataset_->cluster_nodes[trans_edges[i * 2]];
		CNode* node_two = dataset_->cluster_nodes[trans_edges[i * 2 + 1]];

		float node_one_center_x = node_one->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_one_center_y = node_one->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
		float node_two_center_x = node_two->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_two_center_y = node_two->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

		float dir_vec[2], orth_vec[2];
		dir_vec[0] = node_two->center_pos[0] - node_one->center_pos[0];
		dir_vec[1] = node_two->center_pos[1] - node_one->center_pos[1];
		float dir_length = sqrt(pow(dir_vec[0], 2) + pow(dir_vec[1], 2));
		dir_vec[0] /= dir_length;
		dir_vec[1] /= dir_length;
		orth_vec[0] = -1 * dir_vec[1];
		orth_vec[1] = dir_vec[0];

		float begin_pos[2], end_pos[2], arrow_end_pos[2];
		begin_pos[0] = node_one_center_x + node_radius_ * dir_vec[0] * 1.2;
		begin_pos[1] = node_one_center_y + node_radius_ * dir_vec[1] * 1.2;
		end_pos[0] = node_two_center_x - node_radius_ * dir_vec[0] * 1.5;
        end_pos[1] = node_two_center_y - node_radius_ * dir_vec[1] * 1.5;
        arrow_end_pos[0] = node_two_center_x - node_radius_ * dir_vec[0] * 1.2;
        arrow_end_pos[1] = node_two_center_y - node_radius_ * dir_vec[1] * 1.2;

		float value_dis = 0.05;

		cell_ids[0] = trans_points->InsertNextPoint(begin_pos[0] + node_radius_ * orth_vec[0] * value_dis, begin_pos[1] + node_radius_ * orth_vec[1] * value_dis, 0.001);
		cell_ids[1] = trans_points->InsertNextPoint(end_pos[0] + node_radius_ * orth_vec[0] * value_dis, end_pos[1] + node_radius_ * orth_vec[1] * value_dis, 0.001);
		cell_ids[2] = trans_points->InsertNextPoint(end_pos[0] - node_radius_ * orth_vec[0] * value_dis, end_pos[1] - node_radius_ * orth_vec[1] * value_dis, 0.001);
		cell_ids[3] = trans_points->InsertNextPoint(begin_pos[0] - node_radius_ * orth_vec[0] * value_dis, begin_pos[1] - node_radius_ * orth_vec[1] * value_dis, 0.001);
		trans_glyph_polys->InsertNextCell(VTK_POLYGON, 4, cell_ids);

        value_dis = 0.2;
        cell_ids[0] = trans_points->InsertNextPoint(end_pos[0] + node_radius_ * orth_vec[0] * value_dis, end_pos[1] + node_radius_ * orth_vec[1] * value_dis, 0.001);
        cell_ids[1] = trans_points->InsertNextPoint(arrow_end_pos[0], arrow_end_pos[1], 0.001);
		cell_ids[2] = trans_points->InsertNextPoint(end_pos[0] - node_radius_ * orth_vec[0] * value_dis, end_pos[1] - node_radius_ * orth_vec[1] * value_dis, 0.001);
        trans_glyph_polys->InsertNextCell(VTK_POLYGON, 3, cell_ids);
	}

	this->trans_glyph_actors->Modified();
}

void TransMap::UpdateIconActor() {
    if (variable_color_text_.size() < scatter_data_->var_names.size()) {
        for (int i = variable_color_text_.size(); i < scatter_data_->var_names.size(); ++i) {
            vtkTextActor3D* temp_text = vtkTextActor3D::New();
            temp_text->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
	        temp_text->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
	        temp_text->GetTextProperty()->SetBold(true);
	        temp_text->GetTextProperty()->SetFontFamilyToArial();

            QString name = scatter_data_->var_names[i].left(10);
	        temp_text->SetInput(name.toLocal8Bit().data());
	        temp_text->SetPosition(-0.3, -1.6 - (i + 1) * 0.6, 0);
	        temp_text->GetTextProperty()->SetFontSize(20);
	        float scale = 0.2 / 20;
	        temp_text->SetScale(scale, scale, scale);
	        temp_text->Modified();

            variable_color_text_.push_back(temp_text);
            this->glyph_indicator_renderer_->AddActor(temp_text);
        }

        this->glyph_indicator_renderer_->ResetCamera();
        /*double eye_pos[3];
        eye_pos[0] = 0; 
        eye_pos[1] = 0;
        eye_pos[2] = 10;
        float dis = this->glyph_indicator_renderer_->GetActiveCamera()->GetDistance();
        this->glyph_indicator_renderer_->GetActiveCamera()->SetDistance(dis * 5);*/
        this->glyph_indicator_renderer_->GetActiveCamera()->Zoom(1.7);
        this->glyph_indicator_renderer_->GetActiveCamera()->Modified();
    }
    else {
        int item_size = variable_color_text_.size();
        for (int i = scatter_data_->var_names.size(); i < item_size; ++i) {
            this->glyph_indicator_renderer_->RemoveActor(variable_color_text_[i]);
            variable_color_text_[i]->Delete();
        }
        variable_color_text_.resize(scatter_data_->var_names.size());
    }

	var_icon_poly_->Initialize();

    vtkPoints* points = vtkPoints::New();
    vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
	colors->SetNumberOfComponents(4);

	var_icon_poly_->SetPoints(points);
	var_icon_poly_->GetPointData()->SetScalars(colors);

	vtkCellArray* poly_array = vtkCellArray::New();
	var_icon_poly_->SetPolys(poly_array);

    vtkIdType cell_ids[4];
    for (int i = 0; i < scatter_data_->var_names.size(); ++i) {
        float x1 = -1, x2 = -0.45;
        float y1 = -1.6 - i * 0.6 - 0.3;
        float y2 = -1.6 - i * 0.6 - 0.6;

        cell_ids[0] = points->InsertNextPoint(x1, y1, 0.0);
        cell_ids[1] = points->InsertNextPoint(x1, y2, 0.0);
        cell_ids[2] = points->InsertNextPoint(x2, y2, 0.0);
        cell_ids[3] = points->InsertNextPoint(x2, y1, 0.0);
        colors->InsertNextTuple4(scatter_data_->var_colors[i].red(), scatter_data_->var_colors[i].green(), scatter_data_->var_colors[i].blue(), 255);
        colors->InsertNextTuple4(scatter_data_->var_colors[i].red(), scatter_data_->var_colors[i].green(), scatter_data_->var_colors[i].blue(), 255);
        colors->InsertNextTuple4(scatter_data_->var_colors[i].red(), scatter_data_->var_colors[i].green(), scatter_data_->var_colors[i].blue(), 255);
        colors->InsertNextTuple4(scatter_data_->var_colors[i].red(), scatter_data_->var_colors[i].green(), scatter_data_->var_colors[i].blue(), 255);
        
        var_icon_poly_->InsertNextCell(VTK_POLYGON, 4, cell_ids);
    }

    var_icon_actor_->Modified();
}

void TransMap::UpdateDensityActor(std::vector< QColor >& node_colors) {
    this->density_poly_->Initialize();

	vtkPoints* density_points = vtkPoints::New();
	density_poly_->SetPoints(density_points);

	vtkCellArray* density_poly_array = vtkCellArray::New();
	density_poly_->SetPolys(density_poly_array);

    vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
	colors->SetNumberOfComponents(4);
    density_poly_->GetPointData()->SetScalars(colors);

    // 0.1, 0.03
    float radius = 0.03 * scatter_data_->max_pos_range;
    vtkIdType ids[21];
    for (int i = 0; i < node_colors.size(); ++i) {
        float x = scatter_data_->original_point_pos[i][0];
        float y = scatter_data_->original_point_pos[i][1];

        ids[0] = density_points->InsertNextPoint(x, y, -0.0001);
        colors->InsertNextTuple4(node_colors[i].red(), node_colors[i].green(), node_colors[i].blue(), 80);
        
        for (int j = 0; j < 20; ++j) {
            float end_arc = j * 2 * 3.14159 / 19;
            ids[j + 1] = density_points->InsertNextPoint(x + radius * cos(end_arc), y + radius * sin(end_arc), -0.0001);
            colors->InsertNextTuple4(node_colors[i].red(), node_colors[i].green(), node_colors[i].blue(), 10);
        }

        vtkIdType poly_ids[3];
        for (int j = 0; j < 19; ++j) {
            poly_ids[0] = ids[0];
            poly_ids[1] = ids[j + 1];
            poly_ids[2] = ids[j + 2];

            density_poly_->InsertNextCell(VTK_TRIANGLE, 3, poly_ids);
        }
    }

    density_actor_->Modified();
}

void TransMap::UpdateHightlightActor() {
	this->highlight_poly->Initialize();
	int seq_size = seqence_text_actors.size();
	if (seq_size < highlight_node_sequence.size()) {
		for (int i = 0; i < highlight_node_sequence.size() - seq_size; ++i) {
			vtkTextActor3D* actor = vtkTextActor3D::New();
			actor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
			actor->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
			actor->GetTextProperty()->SetBold(true);
			actor->GetTextProperty()->SetFontFamilyToArial();
			this->DefaultRenderer->AddActor(actor);
			seqence_text_actors.push_back(actor);
		}
	}
	else {
		for (int i = seq_size; i >= highlight_node_sequence.size() + 1; --i) {
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
	float x = node_radius_;
	float y = 0;
	vtkIdType cell_ids[5];

	std::list< CNode* >::iterator iter = highlight_node_sequence.begin();
	int hightlight_node_index = 0;
	while (iter != highlight_node_sequence.end()) {
		CNode* temp_node = *iter;
		float node_center_x = temp_node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_center_y = temp_node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

		if (IsLevelOneNode(temp_node)) {
			vtkIdType pre_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.003);
			for (int k = 1; k <= 20; ++k) {
				end_arc = (float)k / 20 * 3.14159 * 2;
				x = node_radius_ * cos(end_arc);
				y = node_radius_ * sin(end_arc);
				vtkIdType current_id = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.003);

				cell_ids[0] = pre_id;
				cell_ids[1] = current_id;
				this->highlight_poly->InsertNextCell(VTK_LINE, 2, cell_ids);

				pre_id = current_id;
			}

			for (int k = 0; k < 21; ++k) color_array->InsertNextTuple3(255, 0.55 * 255, 0.24 * 255);

			char buffer[20];
			itoa(hightlight_node_index, buffer, 10);
			seqence_text_actors[hightlight_node_index]->SetInput(buffer);
			seqence_text_actors[hightlight_node_index]->SetPosition(node_center_x - node_radius_, node_center_y + node_radius_ * 0.5, 0.003);
			seqence_text_actors[hightlight_node_index]->GetTextProperty()->SetFontSize(20);
			float scale = node_radius_ * 0.5 / 30;
			seqence_text_actors[hightlight_node_index]->SetScale(scale, scale, scale);
			seqence_text_actors[hightlight_node_index]->Modified();
			hightlight_node_index++;
		} else {
			vtkIdType pre_id = points->InsertNextPoint(node_center_x + x * 0.4, node_center_y + y * 0.4, 0.003);
			for (int k = 1; k <= 20; ++k) {
				end_arc = (float)k / 20 * 3.14159 * 2;
				x = node_radius_ * cos(end_arc);
				y = node_radius_ * sin(end_arc);
				vtkIdType current_id = points->InsertNextPoint(node_center_x + x * 0.4, node_center_y + y * 0.4, 0.003);

				cell_ids[0] = pre_id;
				cell_ids[1] = current_id;
				this->highlight_poly->InsertNextCell(VTK_LINE, 2, cell_ids);

				pre_id = current_id;
			}

			for (int k = 0; k < 21; ++k) color_array->InsertNextTuple3(255, 0.55 * 255, 0.24 * 255);

			float temp_radius = node_radius_ * 0.2;
			char buffer[20];
			itoa(hightlight_node_index, buffer, 10);
			seqence_text_actors[hightlight_node_index]->SetInput(buffer);
			seqence_text_actors[hightlight_node_index]->SetPosition(node_center_x - temp_radius, node_center_y + temp_radius * 0.5, 0.003);
			seqence_text_actors[hightlight_node_index]->GetTextProperty()->SetFontSize(20);
			float scale = node_radius_ * 0.5 / 30;
			seqence_text_actors[hightlight_node_index]->SetScale(scale, scale, scale);
			seqence_text_actors[hightlight_node_index]->Modified();
			hightlight_node_index++;
		}
		iter++;
	}

	this->highlight_actor->GetProperty()->SetLineStipplePattern(0xF00F);
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

		for (int i = 0; i < node_glyph_actors.size(); ++i) {
			this->DefaultRenderer->AddActor(node_glyph_actors[i]);
			this->node_picker->AddPickList(node_glyph_actors[i]);
		}

		for (int i = 0; i < node_glyph_actors.size(); ++i) {
			this->DefaultRenderer->AddActor(node_glyph_actors[i]);
			this->node_picker->AddPickList(node_glyph_actors[i]);
		}

		this->DefaultRenderer->AddActor(trans_glyph_actors);
        this->DefaultRenderer->AddActor(linked_glyph_actors);
        this->DefaultRenderer->AddActor(density_actor_);
		this->DefaultRenderer->AddActor(this->highlight_actor);
		this->DefaultRenderer->AddActor(this->selection_brush_actor);

        if (this->glyph_indicator_renderer_ != NULL) {
            this->glyph_indicator_renderer_->AddActor(indicator_actor_);
            this->glyph_indicator_renderer_->AddActor(indicator_text_);
            this->glyph_indicator_renderer_->AddActor(var_icon_actor_);

            for (int i = 0; i < variable_color_text_.size(); ++i)
                this->glyph_indicator_renderer_->AddActor(variable_color_text_[i]);
        }

		this->InvokeEvent(vtkCommand::EnableEvent, NULL);
	} else {
		if (!this->Enabled) return;

		this->Enabled = 0;

		this->Interactor->RemoveObserver(this->EventCallbackCommand);

		for (int i = 0; i < node_glyph_actors.size(); ++i) {
			this->DefaultRenderer->RemoveActor(node_glyph_actors[i]);
			this->node_picker->DeletePickList(node_glyph_actors[i]);
		}

		for (int i = 0; i < node_glyph_actors.size(); ++i) {
			this->DefaultRenderer->RemoveActor(node_glyph_actors[i]);
			this->node_picker->DeletePickList(node_glyph_actors[i]);
		}

		this->DefaultRenderer->RemoveActor(trans_glyph_actors);
        this->DefaultRenderer->RemoveActor(linked_glyph_actors);
        this->DefaultRenderer->RemoveActor(density_actor_);
		this->DefaultRenderer->RemoveActor(this->highlight_actor);
		this->DefaultRenderer->RemoveActor(this->selection_brush_actor);

        if (this->glyph_indicator_renderer_ != NULL) {
            this->glyph_indicator_renderer_->RemoveActor(indicator_actor_);
            this->glyph_indicator_renderer_->RemoveActor(indicator_text_);
            this->glyph_indicator_renderer_->RemoveActor(var_icon_actor_);

            for (int i = 0; i < variable_color_text_.size(); ++i)
                this->glyph_indicator_renderer_->RemoveActor(variable_color_text_[i]);

            this->glyph_indicator_renderer_->ResetCamera();
        }

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
	this->UpdateNodeActors();
	this->UpdateHightlightActor();
	this->UpdateTransEdgeActor();
	this->UpdateIconActor();
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
	if (!this->DefaultRenderer) return;

	int x = this->Interactor->GetEventPosition()[0];
	int y = this->Interactor->GetEventPosition()[1];

	switch (this->state)
	{
	case NORMAL:
		break;
	case SELECT_SINGLE_CLUSTER: 
	case SELECT_MULTI_CLUSTERS:
	case SELECT_MINIMUM_PATH: {
		CNode* selected_node = GetSelectedNode(x, y);
		if (selected_node != NULL) OnNodeSelected(selected_node);
		break;
	}
	case SELECT_BRUSHED_PATH_SEQUENCE: {
		this->selection_brush_poly->Initialize();
		vtkPoints* points = vtkPoints::New();
		this->selection_brush_poly->SetPoints(points);
		vtkCellArray* line_array = vtkCellArray::New();
		this->selection_brush_poly->SetLines(line_array);
		this->selection_brush_actor->Modified();

		this->parent_view->update();
		break;
	}
	case SELECT_BRUSHED_CLUSTERS: {
		this->selection_brush_poly->Initialize();
		vtkPoints* points = vtkPoints::New();
		this->selection_brush_poly->SetPoints(points);
		vtkCellArray* line_array = vtkCellArray::New();
		this->selection_brush_poly->SetLines(line_array);
		this->selection_brush_actor->Modified();

		this->parent_view->update();
		break;
	}
	default:
		break;
	}
}

void TransMap::OnRightButtonDown() {
	if (!this->DefaultRenderer) return;

	int x = this->Interactor->GetEventPosition()[0];
	int y = this->Interactor->GetEventPosition()[1];

	CNode* selected_node = GetSelectedNode(x, y);
	if (selected_node != NULL && IsLevelOneNode(selected_node)) {
		double woldpos[4];
		vtkInteractorObserver::ComputeDisplayToWorld(x, y, 0, woldpos);
		float node_center_x = selected_node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_center_y = selected_node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
		float xt = woldpos[0] - node_center_x;
		float yt = woldpos[1] - node_center_y;
		float length = sqrt(pow(xt, 2) + pow(yt, 2));
		float degree = acos(xt / length);
		if (yt < 0) degree = 2 * 3.14159 - degree;
		int temp_index = (int)(degree / (2 * 3.1416) * dataset_->var_num);

        int var_index = axis_order_[temp_index];
        int hindex = -1;
        for (int i = 0; i < this->focus_var_index_.size(); ++i) {
            if (this->focus_var_index_[i] == var_index) {
                hindex = i;
                break;
            }
        }
        if (hindex != -1) {
            for (int i = hindex; i < this->focus_var_index_.size() - 1; ++i)
                this->focus_var_index_[i] = this->focus_var_index_[i + 1];
            this->focus_var_index_.resize(this->focus_var_index_.size() - 1);
        }
        else {
            this->focus_var_index_.push_back(var_index);
        }
	} 

	if (this->focus_var_index_.size() == 0) {
		is_focus_var_fixed_ = false;
		current_node_ = NULL;
	} else {
		is_focus_var_fixed_ = true;
		current_node_ = selected_node;
	}

	this->UpdateNodeActors();
	this->parent_view->update();
}

void TransMap::OnMouseMove() {
	if (!this->DefaultRenderer) return;

	int x = this->Interactor->GetEventPosition()[0];
	int y = this->Interactor->GetEventPosition()[1];
	current_node_ = GetSelectedNode(x, y);

    float pre_highlight = current_highlight_var_index_;

    if (current_node_ == NULL && this->current_highlight_var_index_ != -1) {
        this->UpdateNodeActors();
        this->parent_view->update();
        this->current_highlight_var_index_ = -1;
    }

	current_highlight_var_index_ = -1;
	if (current_node_ != NULL && IsLevelOneNode(current_node_)) {
		double woldpos[4];
		vtkInteractorObserver::ComputeDisplayToWorld(x, y, 0, woldpos);
		float node_center_x = current_node_->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
		float node_center_y = current_node_->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
		float xt = woldpos[0] - node_center_x;
		float yt = woldpos[1] - node_center_y;
		float length = sqrt(pow(xt, 2) + pow(yt, 2));
		float degree = acos(xt / length);
		if (yt < 0) degree = 2 * 3.14159 - degree;
		int temp_index = (int)(degree / (2 * 3.1416) * dataset_->var_num);
		current_highlight_var_index_ = axis_order_[temp_index];
	}
	if (current_highlight_var_index_ != -1) {
		this->UpdateNodeActors();
		this->parent_view->update();

		if (current_highlight_var_index_ != -1) {
			std::vector< float > ranges = dataset_->dataset->original_value_ranges[current_highlight_var_index_];
			float average_value = current_node_->average_values[current_highlight_var_index_]
				* (ranges[1] - ranges[0]) + ranges[0];
			float variance_value = current_node_->variable_variances[current_highlight_var_index_] * (ranges[1] - ranges[0]);

			this->parent_view->ShowTooltip(current_node_->point_count, dataset_->dataset->var_names[current_highlight_var_index_], average_value, variance_value);
		}
		else
			this->parent_view->HideTooltip();
	}
    if (current_highlight_var_index_ != -1) {
        if (!is_hovering_) {
            is_hovering_ = true;
            this->UpdateNodeActors();
            this->parent_view->update();
        }
    }
    else {
        if (is_hovering_) {
            is_hovering_ = false;
            this->UpdateNodeActors();
            this->parent_view->update();
        }
    }

    if (pre_highlight != current_highlight_var_index_) 
        this->parent_view->SetHighlightVarIndex(current_highlight_var_index_);
}

void TransMap::OnMouseMove(int x, int y) {
	if (this->state == WidgetState::SELECT_BRUSHED_CLUSTERS) {
		double display_pos[4];
		display_pos[0] = x;
		display_pos[1] = y;
		display_pos[2] = 0;
		display_pos[3] = 1;
		this->DefaultRenderer->SetDisplayPoint(display_pos);
		this->DefaultRenderer->DisplayToWorld();
		double* world_pos = this->DefaultRenderer->GetWorldPoint();

		this->selection_brush_poly->GetPoints()->InsertNextPoint(world_pos[0], world_pos[1], 0.003);

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
	} else if (this->state == WidgetState::SELECT_BRUSHED_PATH_SEQUENCE) {
		double display_pos[4];
		display_pos[0] = x;
		display_pos[1] = y;
		display_pos[2] = 0;
		display_pos[3] = 1;
		this->DefaultRenderer->SetDisplayPoint(display_pos);
		this->DefaultRenderer->DisplayToWorld();
		double* world_pos = this->DefaultRenderer->GetWorldPoint();

		this->selection_brush_poly->GetPoints()->InsertNextPoint(world_pos[0], world_pos[1], 0.003);

		vtkCellArray* line_array = this->selection_brush_poly->GetLines();
		line_array->Initialize();
		int pnum = this->selection_brush_poly->GetPoints()->GetNumberOfPoints();
		vtkIdType cell_ids[2];
		for (int i = 0; i < pnum - 1; ++i){
			cell_ids[0] = i;
			cell_ids[1] = i + 1;
			this->selection_brush_poly->InsertNextCell(VTK_LINE, 2, cell_ids);
		}

		this->selection_brush_poly->Modified();
		this->selection_brush_actor->Modified();

		this->parent_view->update();
	} 
}

void TransMap::OnMouseReleased(bool is_left_button) {
	if (is_left_button)
		this->OnLeftButtonUp();
	else
		this->OnRightButtonUp();
}

void TransMap::OnLeftButtonUp() {
	if (this->state == WidgetState::SELECT_BRUSHED_CLUSTERS) {
		while (this->highlight_node_sequence.size() > 0) {
			this->highlight_node_sequence.front()->is_highlighted = false;
			this->highlight_node_sequence.pop_front();
		}

		for (int i = 0; i < dataset_->cluster_nodes.size(); ++i) {
			CNode* node = dataset_->cluster_nodes[i];
			float node_center_x = node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
			float node_center_y = node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];
			if (IsInsideSelection(node_center_x, node_center_y)) {
				highlight_node_sequence.push_back(dataset_->cluster_nodes[i]);
				dataset_->cluster_nodes[i]->is_highlighted = true;
			}
		}

		this->UpdateHightlightActor();

		this->selection_brush_poly->Initialize();
		this->selection_brush_actor->Modified();
	} 

	if (this->state == WidgetState::SELECT_BRUSHED_PATH_SEQUENCE) {
		while (this->highlight_node_sequence.size() > 0) {
			this->highlight_node_sequence.front()->is_highlighted = false;
			this->highlight_node_sequence.pop_front();
		}

		std::vector< int > min_dis_index_seq;
		std::vector< int > node_index;

		for (int i = 0; i < dataset_->cluster_nodes.size(); ++i) {
			CNode* node = dataset_->cluster_nodes[i];
			float node_center_x = node->center_pos[0] * (scatter_data_->original_pos_ranges[0][1] - scatter_data_->original_pos_ranges[0][0]) + scatter_data_->original_pos_ranges[0][0];
			float node_center_y = node->center_pos[1] * (scatter_data_->original_pos_ranges[1][1] - scatter_data_->original_pos_ranges[1][0]) + scatter_data_->original_pos_ranges[1][0];

			int min_dis_index = -1;
			float temp_min_dis = 1e10;
			float temp_radius = 0;
			if (node->point_count >= dataset_->min_point_num) {
				temp_radius = node_radius_;
			} else {
				temp_radius = node_radius_ * 0.2;
			}
			
			int count = this->selection_brush_poly->GetPoints()->GetNumberOfPoints();
			double p1[2];
			for (int j = 0; j < count; ++j) {
				double* temp_pos = this->selection_brush_poly->GetPoint(j);
				p1[0] = temp_pos[0] - node_center_x;
				p1[1] = temp_pos[1] - node_center_y;
				double dis = sqrt(pow(p1[0], 2) + pow(p1[1], 2));
				if (dis < temp_radius && dis < temp_min_dis) {
					temp_min_dis = dis;
					min_dis_index = j;
				}
			}
			if (min_dis_index != -1) {
				min_dis_index_seq.push_back(min_dis_index);
				node_index.push_back(i);
			}
		}

		Sort(min_dis_index_seq, node_index);
		for (int i = 0; i < node_index.size(); ++i) {
			dataset_->cluster_nodes[node_index[i]]->is_highlighted = true;
			this->highlight_node_sequence.push_back(dataset_->cluster_nodes[node_index[i]]);
		}

		this->UpdateHightlightActor();

		this->selection_brush_poly->Initialize();
		this->selection_brush_actor->Modified();

		this->GenerateTransEdgeFromHighlight();
		this->UpdateTransEdgeActor();
	}
}

void TransMap::OnNodeSelected(CNode* node) {
	std::list< CNode* >::iterator iter = this->highlight_node_sequence.begin();
	while (iter != this->highlight_node_sequence.end() && (*iter)->id() != node->id()) iter++;
	if (iter != this->highlight_node_sequence.end()) {
		(*iter)->is_highlighted = false;
		this->highlight_node_sequence.erase(iter);
		this->UpdateHightlightActor();

		if (this->state == WidgetState::SELECT_MULTI_CLUSTERS
			|| this->state == WidgetState::SELECT_MINIMUM_PATH 
			|| this->state == WidgetState::SELECT_BRUSHED_PATH_SEQUENCE) {
			this->GenerateTransEdgeFromHighlight();
			this->UpdateTransEdgeActor();
		}

		return;
	}

	switch (this->state)
	{
	case SELECT_SINGLE_CLUSTER: {
		node->is_highlighted = true;
		
		while (this->highlight_node_sequence.size() > 0) {
			this->highlight_node_sequence.front()->is_highlighted = false;
			this->highlight_node_sequence.pop_front();
		}
		this->highlight_node_sequence.push_back(node);
		break;
	}
	case SELECT_MULTI_CLUSTERS: {
		node->is_highlighted = true;
		this->highlight_node_sequence.push_back(node);
		this->GenerateTransEdgeFromHighlight();
		this->UpdateTransEdgeActor();

		break;
	}
	case SELECT_MINIMUM_PATH: {
		if (!IsLevelOneNode(node)) return;

		if (this->highlight_node_sequence.size() >= 2) {
			while (this->highlight_node_sequence.size() > 0) {
				this->highlight_node_sequence.front()->is_highlighted = false;
				this->highlight_node_sequence.pop_front();
			}
		}
		this->highlight_node_sequence.push_back(node);
		if (this->highlight_node_sequence.size() == 2) {
			CNode* first_node = this->highlight_node_sequence.front();
			this->highlight_node_sequence.pop_front();
			CNode* second_node = this->highlight_node_sequence.front();
			this->highlight_node_sequence.pop_front();

			int index_one = this->GetLevelOneNodeIndex(first_node);
			int index_two = this->GetLevelOneNodeIndex(second_node);

			std::vector< int > tour_list;
			this->path_generator_->GenerateMinimumPath(index_one, index_two, tour_list);

			for (int i = 0; i < tour_list.size(); ++i) {
				dataset_->level_one_nodes[tour_list[i]]->is_highlighted = true;
				this->highlight_node_sequence.push_back(dataset_->level_one_nodes[tour_list[i]]);
			}
		}
		this->GenerateTransEdgeFromHighlight();
		this->UpdateTransEdgeActor();
		break;
	}
	default:
		break;
	}

	this->UpdateHightlightActor();
}

void TransMap::GenerateTransEdgeFromHighlight() {
	this->trans_edges.clear();
	if (this->highlight_node_sequence.size() <= 1) return;

	std::list< CNode* >::iterator iter = this->highlight_node_sequence.begin();
	int pre_index = GetClusterNodeIndex(*iter);
	iter++;
	while (iter != this->highlight_node_sequence.end()) {
		int c_index = GetClusterNodeIndex(*iter);
		this->trans_edges.push_back(pre_index);
		this->trans_edges.push_back(c_index);

		pre_index = c_index;
		iter++;
	}
}

int TransMap::GetLevelOneNodeIndex(CNode* node) {
	int node_index = -1;
	for (int i = 0; i < dataset_->level_one_nodes.size(); ++i)
		if (dataset_->level_one_nodes[i] == node) {
			node_index = i;
			break;
		}
	return node_index;
}

int TransMap::GetClusterNodeIndex(CNode* node) {
	int node_index = -1;
	for (int i = 0; i < dataset_->cluster_nodes.size(); ++i)
		if (dataset_->cluster_nodes[i] == node) {
			node_index = i;
			break;
		}
	return node_index;
}

bool TransMap::IsLevelOneNode(CNode* node) {
	int node_index = -1;
	for (int i = 0; i < dataset_->level_one_nodes.size(); ++i)
		if (dataset_->level_one_nodes[i] == node) {
			node_index = i;
			break;
		}
	if (node_index != -1) return true;
	else return false;
}

CNode* TransMap::GetSelectedNode(int x, int y) {
	this->node_picker->Pick(x, y, 1.0, this->DefaultRenderer);
	vtkProp* prop = this->node_picker->GetViewProp();
	if (prop == NULL) return NULL;

	vtkActor* actor = dynamic_cast< vtkActor* >(prop);
	if (actor == NULL) return NULL;

	int node_index = -1;
	for (int i = 0; i < node_glyph_actors.size(); ++i)
		if (node_glyph_actors[i] == actor) {
			node_index = i;
			break;
		}
	if (node_index != -1) 
		return dataset_->cluster_nodes[node_index];

	return NULL;
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

void TransMap::Sort(std::vector< int >& index_one, std::vector< int >& index_two) {
	if (index_one.size() == 0) return;
	for (int i = 0; i < index_one.size() - 1; ++i)
		for (int j = i + 1; j < index_one.size(); ++j)
			if (index_one[i] > index_one[j]){
				float temp_value = index_one[j];
				index_one[j] = index_one[i];
				index_one[i] = temp_value;

				int temp_index = index_two[i];
				index_two[i] = index_two[j];
				index_two[j] = temp_index;
			}
}

void TransMap::SetDensityMapVisibility(bool visibility)
{
    is_densitymap_shown_ = visibility;
    density_actor_->SetVisibility(visibility);

    this->parent_view->update();
}
