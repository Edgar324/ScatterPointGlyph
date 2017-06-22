/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "mean_std_widget.h"

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include "glyph_object.h"
#include "color_mapping_generator.h"

vtkStandardNewMacro(MeanStdWidget);

MeanStdWidget::MeanStdWidget() 
    : GlyphWidget() {
    this->glyph_actor_ = vtkActor::New();
	this->glyph_poly_ = vtkPolyData::New();
	this->glyph_data_mapper_ = vtkPolyDataMapper::New();
	this->glyph_data_mapper_->SetInputData(this->glyph_poly_);
	this->glyph_actor_->SetMapper(this->glyph_data_mapper_);
    this->glyph_actor_->GetProperty()->SetLineWidth(2.0);

    this->highlight_actor_ = vtkActor::New();
	this->highlight_poly_ = vtkPolyData::New();
	this->highlight_data_mapper_ = vtkPolyDataMapper::New();
	this->highlight_data_mapper_->SetInputData(this->highlight_poly_);
	this->highlight_actor_->SetMapper(this->highlight_data_mapper_);
    this->highlight_actor_->GetProperty()->SetLineWidth(3.0);

    this->picker_ = vtkCellPicker::New();
    this->picker_->PickFromListOn();
    this->picker_->AddPickList(this->glyph_actor_);
}

MeanStdWidget::~MeanStdWidget() {

}

void MeanStdWidget::SetEnabled(int enabling) {
    if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

	if (enabling) {
		if (this->Enabled) return;

        if (!this->CurrentRenderer) {
            this->SetCurrentRenderer(this->Interactor->FindPokedRenderer(
                this->Interactor->GetLastEventPosition()[0],
                this->Interactor->GetLastEventPosition()[1]));
            if (this->CurrentRenderer == NULL) return;
        }

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

		this->CurrentRenderer->AddActor(this->glyph_actor_);
        this->CurrentRenderer->AddActor(this->highlight_actor_);

        this->InvokeEvent(vtkCommand::EnableEvent,NULL);
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

        this->Interactor->RemoveObserver(this->EventCallbackCommand);

		this->CurrentRenderer->RemoveActor(this->glyph_actor_);
        this->CurrentRenderer->RemoveActor(this->highlight_actor_);

        this->InvokeEvent(vtkCommand::DisableEvent,NULL);
        this->SetCurrentRenderer(NULL);
	}

	//this->Interactor->Render();
}

void MeanStdWidget::SetRenderingMode(int mode) {
    this->rendering_mode_ = mode;
    this->BuildRepresentation();
}

void MeanStdWidget::BuildRepresentation() {
    if (glyph_object_ == NULL) return;

    if (glyph_object_->point_count() > 5) {
        if (this->rendering_mode_ == 0)
            this->BuildLargeVarGlyph();
        else
            this->BuildLargeMeanGlyph();
    }
    else
        this->BuildSmallGlyph();

    this->glyph_actor_->Modified();
    this->highlight_actor_->Modified();
}


void MeanStdWidget::BuildLargeVarGlyph() {
    vector<QString>& names_ = glyph_object_->names();
    vector<QColor>& colors_ = glyph_object_->colors();
    vector<float>& means_ = glyph_object_->means();
    vector<float>& std_devs_ = glyph_object_->std_devs();
    float node_radius_ = glyph_object_->node_radius();
    float center_x_ = glyph_object_->center_x();
    float center_y_ = glyph_object_->center_y();
    int point_count_ = glyph_object_->point_count();
    int max_point_count_ = glyph_object_->max_point_count();
    float saliency = glyph_object_->saliency();
    float saliency_gray_value = 250 * (1.0 - exp(-3 * (1.0 - saliency)));

    if (names_.size() == 0) return;

    vtkPoints* glyph_points = vtkPoints::New();
	vtkUnsignedCharArray* glyph_colors = vtkUnsignedCharArray::New();
	glyph_colors->SetNumberOfComponents(4);

	glyph_poly_->Initialize();
	glyph_poly_->SetPoints(glyph_points);
	glyph_poly_->GetPointData()->SetScalars(glyph_colors);

	vtkCellArray* glyph_poly_array = vtkCellArray::New();
	glyph_poly_->SetPolys(glyph_poly_array);

	vtkCellArray* glyph_line_array = vtkCellArray::New();
	glyph_poly_->SetLines(glyph_line_array);

	vtkCellArray* glyph_strip_array = vtkCellArray::New();
	glyph_poly_->SetStrips(glyph_strip_array);

    int var_num = names_.size();
    vtkIdType cell_ids[3];

    // paint background
    vector<vtkIdType > background_ids;
	for (int i = 0; i <= SEG_PER_CIRCLE; ++i) {
		float end_arc = i * 3.14159 * 2 / SEG_PER_PIE;
		float x = node_radius_ * cos(end_arc) * 1.0;
		float y = node_radius_ * sin(end_arc) * 1.0;

		background_ids.push_back(glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.001));
				
        glyph_colors->InsertNextTuple4(255, 255, 255, 255);
	}
	glyph_poly_->InsertNextCell(VTK_POLYGON, SEG_PER_CIRCLE + 1, background_ids.data());

    // paint encoding for the point number, 
    /*float point_rate = (float)point_count_ / max_point_count_;
    int seg_num = (SEG_PER_CIRCLE - 1) * point_rate;
    std::vector<vtkIdType > num_inner_ids, num_outer_ids;
    for (int i = 0; i <= seg_num; ++i) {
        float end_arc = -1 * i * 3.14159 * 2 / SEG_PER_CIRCLE + 3.14159 * 0.5;
        float cos_value = cos(end_arc);
        float sin_value = sin(end_arc);
        float x = node_radius_ * cos_value * 1.05;
		float y = node_radius_ * sin_value * 1.05;

        vtkIdType id_one = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.0);
		glyph_colors->InsertNextTuple4(128, 128, 128, 255);
        num_inner_ids.push_back(id_one);

        x = node_radius_ * cos_value * 1.18;
		y = node_radius_ * sin_value * 1.18;
        vtkIdType id_two = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.0);
		glyph_colors->InsertNextTuple4(128, 128, 128, 255);
        num_outer_ids.push_back(id_two);
    }

    
    for (int i = 0; i < num_outer_ids.size() - 1; ++i) {
        cell_ids[0] = num_outer_ids[i];
        cell_ids[1] = num_outer_ids[i + 1];
        cell_ids[2] = num_inner_ids[i];
        glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

        cell_ids[0] = num_outer_ids[i + 1];
        cell_ids[1] = num_inner_ids[i + 1];
        cell_ids[2] = num_inner_ids[i];
        glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
    }*/

    /*vector<vtkIdType> center_cirlce_ids;
	for (int j = 0; j <= SEG_PER_PIE; ++j) {
		float end_arc = j * 3.14159 * 2 / SEG_PER_PIE;
		float x = node_radius_ * 0.05 * cos(end_arc);
		float y = node_radius_ * 0.05 * sin(end_arc);

		center_cirlce_ids.push_back(glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.001));
		glyph_colors->InsertNextTuple4(0, 0, 0, 255);
	}
	glyph_poly_->InsertNextCell(VTK_POLYGON, SEG_PER_PIE + 1, center_cirlce_ids.data());*/

    // paint glyph
    std::vector<vtkIdType > value_pos_ids;
	value_pos_ids.resize(var_num);

	std::vector<int> var_point_ids;
	var_point_ids.resize(var_num * 2);

    std::vector<vtkIdType > mean_value_ids;
    mean_value_ids.resize(var_num + 1);

	for (int i = 0; i < var_num; ++i) {
		float end_arc = i * 3.14159 * 2 / var_num;

		float temp_radius = node_radius_ * (means_[i] - std_devs_[i]) * 0.9 + node_radius_ * 0.1;
		if (temp_radius < node_radius_ * 0.1) temp_radius = node_radius_ * 0.1;
		float x = temp_radius * cos(end_arc);
		float y = temp_radius * sin(end_arc);

		vtkIdType id_one = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
		glyph_colors->InsertNextTuple4(255, 0.55 * 255, 0.24 * 255, 255);

        temp_radius = node_radius_ * means_[i] * 0.9 + node_radius_ * 0.1;
		x = temp_radius * cos(end_arc);
		y = temp_radius * sin(end_arc);

		value_pos_ids[i] = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
		glyph_colors->InsertNextTuple4(255, 0.55 * 255, 0.24 * 255, 255);

        mean_value_ids[i] = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 3);
		glyph_colors->InsertNextTuple4(64, 64, 64, 255);

		temp_radius = node_radius_ * (means_[i] + std_devs_[i]) * 0.9 + node_radius_ * 0.1;
		if (temp_radius > node_radius_) temp_radius = node_radius_;
		x = temp_radius * cos(end_arc);
		y = temp_radius * sin(end_arc);

		vtkIdType id_two = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
		glyph_colors->InsertNextTuple4(255, 0.55 * 255, 0.24 * 255, 255);

		var_point_ids[i * 2] = id_one;
		var_point_ids[i * 2 + 1] = id_two;
	}

    mean_value_ids[mean_value_ids.size() - 1] = mean_value_ids[0];
    glyph_poly_->InsertNextCell(VTK_POLY_LINE, mean_value_ids.size(), mean_value_ids.data());


	for (int i = 0; i < var_num; ++i) {
		cell_ids[0] = var_point_ids[2 * i];
		cell_ids[1] = value_pos_ids[i];
		cell_ids[2] = var_point_ids[2 * ((i + 1) % var_num)];
		glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
        cell_ids[0] = value_pos_ids[(i + 1) % var_num];
		cell_ids[1] = value_pos_ids[i];
		cell_ids[2] = var_point_ids[2 * ((i + 1) % var_num)];
		glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

		cell_ids[0] = var_point_ids[2 * i + 1];
		cell_ids[1] = value_pos_ids[i];
		cell_ids[2] = var_point_ids[2 * ((i + 1) % var_num) + 1];
		glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
        cell_ids[0] = value_pos_ids[(i + 1) % var_num];
		cell_ids[1] = value_pos_ids[i];
		cell_ids[2] = var_point_ids[2 * ((i + 1) % var_num) + 1];
		glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
	}

    

    vtkIdType center_id = glyph_points->InsertNextPoint(center_x_, center_y_, 1.0);
	glyph_colors->InsertNextTuple4(220, 220, 220, 255);
	for (int j = 0; j < var_num; ++j) {
		float end_arc = j * 3.14159 * 2 / var_num;
		float x = node_radius_ * cos(end_arc);
		float y = node_radius_ * sin(end_arc);

		vtkIdType pre_id = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 1.0);
		glyph_colors->InsertNextTuple4(220, 220, 220, 255);

		cell_ids[0] = center_id;
		cell_ids[1] = pre_id;
		glyph_poly_->InsertNextCell(VTK_LINE, 2, cell_ids);
	}

    std::vector<vtkIdType > range_outer_ids, range_inner_ids, saliency_ids;
    for (int j = 0; j < SEG_PER_CIRCLE; ++j) {
        float end_arc = -1 * j * 3.14159 * 2 / (SEG_PER_CIRCLE - 1);
        float cos_value = cos(end_arc);
        float sin_value = sin(end_arc);

        float x = node_radius_ * cos_value;
        float y = node_radius_ * sin_value;
        vtkIdType id_three = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 1.0);
		glyph_colors->InsertNextTuple4(220, 220, 220, 255);
        range_outer_ids.push_back(id_three);
    }

    glyph_poly_->InsertNextCell(VTK_POLY_LINE, range_outer_ids.size(), range_outer_ids.data());

    // update highlight actor
    vtkPoints* highlight_points = vtkPoints::New();
	vtkUnsignedCharArray* highlight_colors = vtkUnsignedCharArray::New();
	highlight_colors->SetNumberOfComponents(4);

	highlight_poly_->Initialize();
	highlight_poly_->SetPoints(highlight_points);
	highlight_poly_->GetPointData()->SetScalars(highlight_colors);

	vtkCellArray* highlight_line_array = vtkCellArray::New();
	highlight_poly_->SetLines(highlight_line_array);

    // Add highlight index
    int highlight_index = glyph_object_->highlight_index();
    if (highlight_index != -1) {
        float highlight_val = glyph_object_->highlight_val();
        float begin_arc = (highlight_index - 0.5) * PIE_VAL * 2 / var_num;
        float end_arc = (highlight_index + 0.5) * PIE_VAL * 2 / var_num;
        float step_arc = (end_arc - begin_arc) / (SEG_PER_PIE - 1);

		vector<int> comp_ids;
		comp_ids.resize(SEG_PER_PIE);

		float temp_arc = begin_arc;
		for (int k = 0; k < SEG_PER_PIE; ++k) {
			float cos_value = cos(temp_arc);
			float sin_value = sin(temp_arc);

			// insert the average value
			float x = node_radius_ * cos_value;
			float y = node_radius_ * sin_value;
			comp_ids[k] = highlight_points->InsertNextPoint(center_x_ + x, center_y_ + y, 3.0);
			highlight_colors->InsertNextTuple4(255, 50, 50, 200);

			temp_arc += step_arc;
		}

        vtkIdType line_ids[2];
		for (int k = 0; k < SEG_PER_PIE - 1; ++k) {
			line_ids[0] = comp_ids[k];
			line_ids[1] = comp_ids[k + 1];
			highlight_poly_->InsertNextCell(VTK_LINE, 2, line_ids);
		}
	}

    // paint expandable indicator
    if (glyph_object_->is_expandable()) {
        vtkIdType line_ids[2];

        float x1 = center_x_ - node_radius_ * 1.4;
        float x2 = center_x_ - node_radius_ * 1.0;
        float x3 = center_x_ - node_radius_ * 1.2;
        float y1 = center_y_ + node_radius_ * 0.8;
        float y2 = center_y_ + node_radius_ * 1.2;
        float y3 = center_y_ + node_radius_ * 1.0;

        line_ids[0] = highlight_points->InsertNextPoint(x1, y3, 0);
        line_ids[1] = highlight_points->InsertNextPoint(x2, y3, 0);
        highlight_colors->InsertNextTuple4(0, 0, 0, 255);
        highlight_colors->InsertNextTuple4(0, 0, 0, 255);
        highlight_poly_->InsertNextCell(VTK_LINE, 2, line_ids);

        line_ids[0] = highlight_points->InsertNextPoint(x3, y1, 0);
        line_ids[1] = highlight_points->InsertNextPoint(x3, y2, 0);
        highlight_colors->InsertNextTuple4(0, 0, 0, 255);
        highlight_colors->InsertNextTuple4(0, 0, 0, 255);
        highlight_poly_->InsertNextCell(VTK_LINE, 2, line_ids);
    }

    // add selection indicator
    if (glyph_object_->is_selected()) {
        float begin_arc = 0;
        float end_arc = PIE_VAL * 2;
        float step_arc = (end_arc - begin_arc) / (SEG_PER_CIRCLE - 1);

	    vector<vtkIdType> comp_ids;
	    comp_ids.resize(SEG_PER_CIRCLE);

	    float temp_arc = begin_arc;
	    for (int k = 0; k < SEG_PER_CIRCLE; ++k) {
		    float cos_value = cos(temp_arc);
		    float sin_value = sin(temp_arc);

		    // insert the average value
		    float temp_radius = node_radius_ * 1.0;
		    float x = temp_radius * cos_value;
		    float y = temp_radius * sin_value;
		    comp_ids[k] = highlight_points->InsertNextPoint(center_x_ + x, center_y_ + y, 3);
            highlight_colors->InsertNextTuple4(255, 0, 0, 255);

		    temp_arc += step_arc;
	    }

        highlight_poly_->InsertNextCell(VTK_POLY_LINE, comp_ids.size(), comp_ids.data());
    }
}


void MeanStdWidget::BuildLargeMeanGlyph() {
    vector<QString>& names_ = glyph_object_->names();
    vector<QColor>& colors_ = glyph_object_->colors();
    vector<float>& means_ = glyph_object_->means();
    vector<float>& std_devs_ = glyph_object_->std_devs();
    vector<float>& bias = glyph_object_->bias();
    vector<float>& dev_ranges = glyph_object_->dev_ranges();
    float node_radius_ = glyph_object_->node_radius();
    float center_radius = 0.5 * node_radius_;
    float half_size = 0.5 * node_radius_;
    float center_x_ = glyph_object_->center_x();
    float center_y_ = glyph_object_->center_y();
    int point_count_ = glyph_object_->point_count();
    int max_point_count_ = glyph_object_->max_point_count();

    if (names_.size() == 0) return;

    vtkPoints* glyph_points = vtkPoints::New();
	vtkUnsignedCharArray* glyph_colors = vtkUnsignedCharArray::New();
	glyph_colors->SetNumberOfComponents(4);

	glyph_poly_->Initialize();
	glyph_poly_->SetPoints(glyph_points);
	glyph_poly_->GetPointData()->SetScalars(glyph_colors);

	vtkCellArray* glyph_poly_array = vtkCellArray::New();
	glyph_poly_->SetPolys(glyph_poly_array);

	vtkCellArray* glyph_line_array = vtkCellArray::New();
	glyph_poly_->SetLines(glyph_line_array);

	vtkCellArray* glyph_strip_array = vtkCellArray::New();
	glyph_poly_->SetStrips(glyph_strip_array);

    int var_num = names_.size();
    vtkIdType cell_ids[3];

    // paint background
    vector<vtkIdType > background_ids;
	for (int j = 0; j <= SEG_PER_CIRCLE; ++j) {
		float end_arc = j * 3.14159 * 2 / SEG_PER_PIE;
		float x = node_radius_ * cos(end_arc) * 1.0;
		float y = node_radius_ * sin(end_arc) * 1.0;

		background_ids.push_back(glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.001));
				
        glyph_colors->InsertNextTuple4(255, 255, 255, 255);
	}
	glyph_poly_->InsertNextCell(VTK_POLYGON, SEG_PER_CIRCLE + 1, background_ids.data());

    // paint encoding for the point number, 
    /*float point_rate = (float)point_count_ / max_point_count_;
    int seg_num = (SEG_PER_CIRCLE - 1) * point_rate;
    std::vector<vtkIdType > num_inner_ids, num_outer_ids;
    for (int i = 0; i <= seg_num; ++i) {
        float end_arc = -1 * i * 3.14159 * 2 / SEG_PER_CIRCLE + 3.14159 * 0.5;
        float cos_value = cos(end_arc);
        float sin_value = sin(end_arc);
        float x = node_radius_ * cos_value * 1.05;
		float y = node_radius_ * sin_value * 1.05;

        vtkIdType id_one = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.0);
		glyph_colors->InsertNextTuple4(128, 128, 128, 255);
        num_inner_ids.push_back(id_one);

        x = node_radius_ * cos_value * 1.18;
		y = node_radius_ * sin_value * 1.18;
        vtkIdType id_two = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.0);
		glyph_colors->InsertNextTuple4(128, 128, 128, 255);
        num_outer_ids.push_back(id_two);
    }

    
    for (int i = 0; i < num_outer_ids.size() - 1; ++i) {
        cell_ids[0] = num_outer_ids[i];
        cell_ids[1] = num_outer_ids[i + 1];
        cell_ids[2] = num_inner_ids[i];
        glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

        cell_ids[0] = num_outer_ids[i + 1];
        cell_ids[1] = num_inner_ids[i + 1];
        cell_ids[2] = num_inner_ids[i];
        glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
    }*/

    vtkIdType center_id = glyph_points->InsertNextPoint(center_x_, center_y_, 1.0);
	glyph_colors->InsertNextTuple4(220, 220, 220, 255);
	for (int j = 0; j < var_num; ++j) {
		float end_arc = j * 3.14159 * 2 / var_num;
		float x = node_radius_ * cos(end_arc);
		float y = node_radius_ * sin(end_arc);

		vtkIdType pre_id = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 1.0);
		glyph_colors->InsertNextTuple4(220, 220, 220, 255);

		cell_ids[0] = center_id;
		cell_ids[1] = pre_id;
		glyph_poly_->InsertNextCell(VTK_LINE, 2, cell_ids);
	}

    std::vector<vtkIdType > range_outer_ids, range_inner_ids, saliency_ids;
    for (int j = 0; j < SEG_PER_CIRCLE; ++j) {
        float end_arc = -1 * j * 3.14159 * 2 / (SEG_PER_CIRCLE - 1);
        float cos_value = cos(end_arc);
        float sin_value = sin(end_arc);

        float x = node_radius_ * cos_value;
        float y = node_radius_ * sin_value;
        vtkIdType id_three = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 1.0);
		glyph_colors->InsertNextTuple4(220, 220, 220, 255);
        range_outer_ids.push_back(id_three);
    }

    glyph_poly_->InsertNextCell(VTK_POLY_LINE, range_outer_ids.size(), range_outer_ids.data());

    // render mean bias glyph
    int begin_index = 0;
    while (begin_index < bias.size() && bias[begin_index] * bias[(begin_index + 1) % var_num] > 0) begin_index++;
    // step to the next point, make sure the begin_index is the point next to the crossing
    begin_index++;

    if (begin_index >= bias.size()) {
        vector<QPointF> region_points;
        for (int i = 0; i < bias.size(); i++) {
            region_points.push_back(QPointF(i, bias[i]));
        }
        region_points.push_back(QPointF(bias.size(), bias[0]));
        
        vector<QPointF> new_region_points;
        MakeCubic(region_points, SEG_PER_CIRCLE, new_region_points);
        region_points = new_region_points;

        vector<QPointF> zero_points;
        float delta = 0.01;
        if (region_points[1].ry() >= 0) delta = -0.01;

        for (int i = 0; i < region_points.size(); i++) {
            zero_points.push_back(QPointF(region_points[i].rx(), delta));
        }
        
        QColor region_color(49,130,189, 255);
        if (region_points[1].ry() >= 0) region_color = QColor(250,159,181, 255);

        vector<vtkIdType> poly_ids(region_points.size() * 2);
        for (int i = 0; i < region_points.size(); i++) {
            float x = (center_radius + region_points[i].ry() * half_size) * cos(2 * PIE_VAL * region_points[i].rx() / var_num);
            float y = (center_radius + region_points[i].ry() * half_size) * sin(2 * PIE_VAL * region_points[i].rx() / var_num);
            poly_ids[2 * i + 1] = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
            glyph_colors->InsertNextTuple4(region_color.red(), region_color.green(), region_color.blue(), 255);

            x = (center_radius + zero_points[i].ry() * half_size) * cos(2 * PIE_VAL * zero_points[i].rx() / var_num);
            y = (center_radius + zero_points[i].ry() * half_size) * sin(2 * PIE_VAL * zero_points[i].rx() / var_num);
            poly_ids[2 * i] = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
            glyph_colors->InsertNextTuple4(region_color.red(), region_color.green(), region_color.blue(), 255);
        }
        glyph_poly_->InsertNextCell(VTK_TRIANGLE_STRIP, (int)poly_ids.size(), poly_ids.data());
    } else {
        vector<float> rotated_bias;
        for (int i = 0; i < bias.size(); i++) 
            rotated_bias.push_back(bias[(i + begin_index) % var_num]);
        rotated_bias.push_back(0);

        int current_index = 0;
        vector<QPointF> region_points;
        QPointF prePoint(-0.5, 0);

        region_points.push_back(prePoint);
        while (current_index < rotated_bias.size()) {
            float currentY = rotated_bias[current_index];
            float currentX = current_index;

            if (prePoint.ry() * currentY < 0 || currentY == 0) {
                vector<QPointF> zero_points;
                float beginX = currentX - 0.5;
                region_points.push_back(QPointF(beginX, 0));

                QColor region_color(49,130,189, 255);
                if (region_points[1].ry() >= 0) region_color = QColor(250,159,181, 255);

                vector<QPointF> new_region_points;
                int point_num = (beginX - region_points[0].rx()) / var_num * SEG_PER_CIRCLE;
                MakeCubic(region_points, point_num, new_region_points);

                region_points = new_region_points;

                float delta = 0.01;
                if (region_points[1].ry() >= 0) delta = -0.01;

                for (int i = 0; i < region_points.size(); i++) {
                    zero_points.push_back(QPointF(region_points[i].rx(), delta));
                }

                for (int i = 0; i < region_points.size(); i++) {
                    region_points[i].setX(region_points[i].rx() + begin_index);
                    zero_points[i].setX(zero_points[i].rx() + begin_index);
                }

                // insert the region to the glyph
                vector<vtkIdType> poly_ids(region_points.size() * 2);
                if (region_points[1].ry() > 0) {
                    for (int i = 0; i < region_points.size(); i++) {
                        float x = (center_radius + region_points[i].ry() * half_size) * cos(2 * PIE_VAL * region_points[i].rx() / var_num);
                        float y = (center_radius + region_points[i].ry() * half_size) * sin(2 * PIE_VAL * region_points[i].rx() / var_num);
                        poly_ids[2 * i + 1] = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
                        glyph_colors->InsertNextTuple4(region_color.red(), region_color.green(), region_color.blue(), 255);

                        x = (center_radius + zero_points[i].ry() * half_size) * cos(2 * PIE_VAL * zero_points[i].rx() / var_num);
                        y = (center_radius + zero_points[i].ry() * half_size) * sin(2 * PIE_VAL * zero_points[i].rx() / var_num);
                        poly_ids[2 * i] = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
                        glyph_colors->InsertNextTuple4(region_color.red(), region_color.green(), region_color.blue(), 255);
                    }
                } else {
                    for (int i = 0; i < region_points.size(); i++) {
                        float x = (center_radius + region_points[i].ry() * half_size) * cos(2 * PIE_VAL * region_points[i].rx() / var_num);
                        float y = (center_radius + region_points[i].ry() * half_size) * sin(2 * PIE_VAL * region_points[i].rx() / var_num);
                        poly_ids[2 * i] = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
                        glyph_colors->InsertNextTuple4(region_color.red(), region_color.green(), region_color.blue(), 255);

                        x = (center_radius + zero_points[i].ry() * half_size) * cos(2 * PIE_VAL * zero_points[i].rx() / var_num);
                        y = (center_radius + zero_points[i].ry() * half_size) * sin(2 * PIE_VAL * zero_points[i].rx() / var_num);
                        poly_ids[2 * i + 1] = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 2);
                        glyph_colors->InsertNextTuple4(region_color.red(), region_color.green(), region_color.blue(), 255);
                    }
                }
                glyph_poly_->InsertNextCell(VTK_TRIANGLE_STRIP, (int)poly_ids.size(), poly_ids.data());

                prePoint.setX(currentX - 0.5);
                prePoint.setY(0);

                region_points.clear();
                region_points.push_back(prePoint);

                if (current_index + 1 == rotated_bias.size()) break;
            } else {
                /*int point_num = SEG_PER_PIE * (currentX - region_points[0].rx());
                this->CubicTo(region_points[region_points.size() - 1], QPointF(currentX, currentY), point_num, region_points);*/
                region_points.push_back(QPointF(currentX, currentY));
                prePoint.setX(currentX);
                prePoint.setY(currentY);
                current_index++;
            }
        }
    }

    glyph_poly_->Modified();

    vtkPoints* highlight_points = vtkPoints::New();
	vtkUnsignedCharArray* highlight_colors = vtkUnsignedCharArray::New();
	highlight_colors->SetNumberOfComponents(4);

	highlight_poly_->Initialize();
	highlight_poly_->SetPoints(highlight_points);
	highlight_poly_->GetPointData()->SetScalars(highlight_colors);

	vtkCellArray* highlight_line_array = vtkCellArray::New();
	highlight_poly_->SetLines(highlight_line_array);

    // Add highlight index
    int highlight_index = glyph_object_->highlight_index();
    if (highlight_index != -1) {
        float highlight_val = glyph_object_->highlight_val();
        float begin_arc = (highlight_index - 0.5) * PIE_VAL * 2 / var_num;
        float end_arc = (highlight_index + 0.5) * PIE_VAL * 2 / var_num;
        float step_arc = (end_arc - begin_arc) / (SEG_PER_PIE - 1);

		vector<int> comp_ids;
		comp_ids.resize(SEG_PER_PIE);

		float temp_arc = begin_arc;
		for (int k = 0; k < SEG_PER_PIE; ++k) {
			float cos_value = cos(temp_arc);
			float sin_value = sin(temp_arc);

			// insert the average value
			float x = node_radius_ * cos_value;
			float y = node_radius_ * sin_value;
			comp_ids[k] = highlight_points->InsertNextPoint(center_x_ + x, center_y_ + y, 3);
			highlight_colors->InsertNextTuple4(255, 50, 50, 200);

			temp_arc += step_arc;
		}

        vtkIdType line_ids[2];
		for (int k = 0; k < SEG_PER_PIE - 1; ++k) {
			line_ids[0] = comp_ids[k];
			line_ids[1] = comp_ids[k + 1];
			highlight_poly_->InsertNextCell(VTK_LINE, 2, line_ids);
		}
	}

    // paint expandable indicator
    if (glyph_object_->is_expandable()) {
        vtkIdType line_ids[2];

        float x1 = center_x_ - node_radius_ * 1.4;
        float x2 = center_x_ - node_radius_ * 1.0;
        float x3 = center_x_ - node_radius_ * 1.2;
        float y1 = center_y_ + node_radius_ * 0.8;
        float y2 = center_y_ + node_radius_ * 1.2;
        float y3 = center_y_ + node_radius_ * 1.0;

        line_ids[0] = highlight_points->InsertNextPoint(x1, y3, 0.0000001);
        line_ids[1] = highlight_points->InsertNextPoint(x2, y3, 0.0000001);
        highlight_colors->InsertNextTuple4(0, 0, 0, 255);
        highlight_colors->InsertNextTuple4(0, 0, 0, 255);
        highlight_poly_->InsertNextCell(VTK_LINE, 2, line_ids);

        line_ids[0] = highlight_points->InsertNextPoint(x3, y1, 0.0000001);
        line_ids[1] = highlight_points->InsertNextPoint(x3, y2, 0.0000001);
        highlight_colors->InsertNextTuple4(0, 0, 0, 255);
        highlight_colors->InsertNextTuple4(0, 0, 0, 255);
        highlight_poly_->InsertNextCell(VTK_LINE, 2, line_ids);
    }

    // add selection indicator
    if (glyph_object_->is_selected()) {
        float begin_arc = 0;
        float end_arc = PIE_VAL * 2;
        float step_arc = (end_arc - begin_arc) / (SEG_PER_CIRCLE - 1);

	    vector<vtkIdType> comp_ids;
	    comp_ids.resize(SEG_PER_CIRCLE);

	    float temp_arc = begin_arc;
	    for (int k = 0; k < SEG_PER_CIRCLE; ++k) {
		    float cos_value = cos(temp_arc);
		    float sin_value = sin(temp_arc);

		    // insert the average value
		    float temp_radius = node_radius_ * 1.0;
		    float x = temp_radius * cos_value;
		    float y = temp_radius * sin_value;
		    comp_ids[k] = highlight_points->InsertNextPoint(center_x_ + x, center_y_ + y, 3);
            highlight_colors->InsertNextTuple4(255, 0, 0, 255);

		    temp_arc += step_arc;
	    }

        highlight_poly_->InsertNextCell(VTK_POLY_LINE, comp_ids.size(), comp_ids.data());
    }

    highlight_poly_->Modified();
}

void MeanStdWidget::BuildSmallGlyph() {
    vtkPoints* glyph_points = vtkPoints::New();
	vtkUnsignedCharArray* glyph_colors = vtkUnsignedCharArray::New();
	glyph_colors->SetNumberOfComponents(4);

	glyph_poly_->Initialize();
	glyph_poly_->SetPoints(glyph_points);
	glyph_poly_->GetPointData()->SetScalars(glyph_colors);

	vtkCellArray* glyph_poly_array = vtkCellArray::New();
	glyph_poly_->SetPolys(glyph_poly_array);

	vtkCellArray* glyph_line_array = vtkCellArray::New();
	glyph_poly_->SetLines(glyph_line_array);

	vtkCellArray* glyph_strip_array = vtkCellArray::New();
	glyph_poly_->SetStrips(glyph_strip_array);

    float node_radius_ = glyph_object_->node_radius();
    float center_x_ = glyph_object_->center_x();
    float center_y_ = glyph_object_->center_y();
	float temp_radius = node_radius_ * 0.2;
	float end_arc = 0;
	float x = temp_radius;
	float y = 0;

	vtkIdType center_id = glyph_points->InsertNextPoint(center_x_, center_y_, 0.001);
	vtkIdType pre_id = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.001);

    vtkIdType cell_ids[3];
	for (int k = 1; k <= 20; ++k) {
		end_arc = (float)k / 20 * 3.14159 * 2;
		x = temp_radius * cos(end_arc);
		y = temp_radius * sin(end_arc);
		vtkIdType current_id = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.001);

		cell_ids[0] = pre_id;
		cell_ids[1] = current_id;
		cell_ids[2] = center_id;
		glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

		pre_id = current_id;
	}

	for (int k = 0; k < 22; ++k)
		glyph_colors->InsertNextTuple4(255, 140, 61, 255);

    vtkPoints* highlight_points = vtkPoints::New();
	vtkUnsignedCharArray* highlight_colors = vtkUnsignedCharArray::New();
	highlight_colors->SetNumberOfComponents(4);

	highlight_poly_->Initialize();
	highlight_poly_->SetPoints(highlight_points);
	highlight_poly_->GetPointData()->SetScalars(highlight_colors);

	vtkCellArray* highlight_line_array = vtkCellArray::New();
	highlight_poly_->SetLines(highlight_line_array);

    // add selection indicator
    if (glyph_object_->is_selected()) {
        float begin_arc = 0;
        float end_arc = PIE_VAL * 2;
        float step_arc = (end_arc - begin_arc) / (SEG_PER_PIE - 1);

	    vector<int> comp_ids;
	    comp_ids.resize(SEG_PER_PIE);

	    float temp_arc = begin_arc;
	    for (int k = 0; k < SEG_PER_PIE; ++k) {
		    float cos_value = cos(temp_arc);
		    float sin_value = sin(temp_arc);

		    // insert the average value
		    float temp_radius = node_radius_ * 0.2;
		    float x = temp_radius * cos_value;
		    float y = temp_radius * sin_value;
		    comp_ids[k] = highlight_points->InsertNextPoint(center_x_ + x, center_y_ + y, 1e-3);
            highlight_colors->InsertNextTuple4(255, 0, 0, 255);

		    temp_arc += step_arc;
	    }

        vtkIdType line_ids[2];
	    for (int k = 0; k < SEG_PER_PIE - 1; ++k) {
		    line_ids[0] = comp_ids[k];
		    line_ids[1] = comp_ids[k + 1];
		    highlight_poly_->InsertNextCell(VTK_LINE, 2, line_ids);
	    }
    }
}

void MeanStdWidget::MakeCubic(vector<QPointF>& points, int point_num, vector<QPointF>& splines) {
    vector<QPointF> control_points;
    vector<float> control_t;
    int curve_level = 3;

    for (int i = 0; i < curve_level - 1; i++) {
        control_t.push_back(points[0].rx());
        control_points.push_back(points[0]);
    }

    for (int i = 0; i < points.size() - 1; i++) {
        control_points.push_back(points[i]);
        control_t.push_back(points[i].rx());

        float x = points[i].rx() + 0.35 * (points[i + 1].rx() - points[i].rx());
        control_points.push_back(QPointF(x, points[i].ry()));
        control_t.push_back(x);

        x = points[i].rx() + 0.75 * (points[i + 1].rx() - points[i].rx());
        control_points.push_back(QPointF(x, points[i + 1].ry()));
        control_t.push_back(x);


        control_points.push_back(points[i + 1]);
        control_t.push_back(points[i + 1].rx());
    }

    for (int i = 0; i < curve_level - 1; i++) {
        control_t.push_back(points[points.size() - 1].rx());
        control_points.push_back(points[points.size() - 1]);
    }

    
    if (points.size() == 3) curve_level = 2;
    if (points.size() == 2) {
        splines.assign(points.begin(), points.end());
        return;
    }

    float end_x = control_t[control_t.size() - 1];
    float begin_x = control_t[0];

    vector<QPointF> current_points, temp_points;
    for (int i = 0; i < point_num; i++) {
        float t_value = (float)i / (point_num - 1) * (end_x - begin_x) + begin_x;
        int range = 0;
        while (range + 1 < control_t.size() && t_value >= control_t[range + 1]) range++;

        current_points = control_points;
        for (int r = 1; r < curve_level; ++r) {
            temp_points = current_points;
            for (int j = range - curve_level + r + 1; j <= range; j++) {
                float scale = (t_value - control_t[j]);
                scale /= (control_t[j + curve_level - r] - control_t[j]);

                current_points[j] = scale * temp_points[j + 1] + (1.0 - scale) * temp_points[j];
            }
        }
        splines.push_back(current_points[range]);
    }
}

void MeanStdWidget::OnMouseMove() {
    int pos_x = this->Interactor->GetEventPosition()[0];
    int pos_y = this->Interactor->GetEventPosition()[1];

    if (!this->CurrentRenderer || !this->CurrentRenderer->IsInViewport(pos_x, pos_y)) {
        this->glyph_state_ = GlyphWidget::NORMAL;
        return;
    }

    this->picker_->Pick(pos_x, pos_y, 1.0, this->CurrentRenderer);
	vtkAssemblyPath* path = this->picker_->GetPath();
	if (path == NULL) return;

    if (this->glyph_state_ == GlyphWidget::NORMAL) {
        double woldpos[4];
	    this->ComputeDisplayToWorld(pos_x, pos_y, 0, woldpos);
	    float xt = woldpos[0] - glyph_object_->center_x();
	    float yt = woldpos[1] - glyph_object_->center_y();
	    float length = sqrt(pow(xt, 2) + pow(yt, 2));
	    float degree = acos(xt / length);
	    if (yt < 0) degree = 2 * 3.14159 - degree;

	    int highlight_index = (int)(degree / (2 * 3.1416) * glyph_object_->var_num() + 0.5);
        if (highlight_index >= glyph_object_->var_num()) highlight_index = 0;
        if (highlight_index != glyph_object_->highlight_index()) {
            glyph_object_->set_highlight(highlight_index, glyph_object_->means()[highlight_index]);
            if (this->rendering_widget_) {
                this->rendering_widget_->HighlightModified(this->glyph_object_);
            }
        }
    } /*else if (this->glyph_state_ == GlyphWidget::SELECTION) {
        this->glyph_object_->set_is_selected(true);
        if (this->rendering_widget_) {
            this->rendering_widget_->SelectionModified(this->glyph_object_);
        }
    }*/

    this->EventCallbackCommand->SetAbortFlag(1);
    this->Interactor->Render();
} 

void MeanStdWidget::OnLeftButtonDown() {
    int pos_x = this->Interactor->GetEventPosition()[0];
    int pos_y = this->Interactor->GetEventPosition()[1];

    if (!this->CurrentRenderer || !this->CurrentRenderer->IsInViewport(pos_x, pos_y)) {
        this->glyph_state_ = GlyphWidget::NORMAL;
        return;
    }

    //this->glyph_state_ = GlyphWidget::SELECTION;
    glyph_object_->set_highlight(-1, 0);
    if (this->rendering_widget_) {
        this->rendering_widget_->HighlightModified(this->glyph_object_);
    }

    this->picker_->Pick(pos_x, pos_y, 1.0, this->CurrentRenderer);
	vtkAssemblyPath* path = this->picker_->GetPath();
	if (path == NULL) return;

    this->glyph_object_->set_is_selected(!this->glyph_object_->is_selected());
    if (this->rendering_widget_) {
        this->rendering_widget_->SelectionModified(this->glyph_object_);
    }
}

void MeanStdWidget::OnLeftButtonUp() {
    this->glyph_state_ = GlyphWidget::NORMAL;
}

void MeanStdWidget::OnRightButtonDown() {

}

void MeanStdWidget::OnRightButtonUp() {

}