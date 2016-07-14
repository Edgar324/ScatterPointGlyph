/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "radar_band_widget_no_saliency.h"
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
#include "glyph_object.h"

vtkStandardNewMacro(RadarBandWidgetNoSaliency);

RadarBandWidgetNoSaliency::RadarBandWidgetNoSaliency() 
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

RadarBandWidgetNoSaliency::~RadarBandWidgetNoSaliency() {

}

void RadarBandWidgetNoSaliency::SetEnabled(int enabling) {
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

void RadarBandWidgetNoSaliency::BuildRepresentation() {
    if (glyph_object_ == NULL) return;

    if (glyph_object_->point_count() > 5) 
        this->BuildLargeGlyph();
    else
        this->BuildSmallGlyph();

    this->glyph_actor_->Modified();
    this->highlight_actor_->Modified();
}

void RadarBandWidgetNoSaliency::BuildLargeGlyph() {
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

    // paint background
    vector<vtkIdType > background_ids;
	for (int j = 0; j <= SEG_PER_CIRCLE; ++j) {
		float end_arc = j * 3.14159 * 2 / SEG_PER_PIE;
		float x = node_radius_ * cos(end_arc) * 1.0;
		float y = node_radius_ * sin(end_arc) * 1.0;

		background_ids.push_back(glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.000));
				
        glyph_colors->InsertNextTuple4(255, 255, 255, 255);
	}
	glyph_poly_->InsertNextCell(VTK_POLYGON, SEG_PER_CIRCLE + 1, background_ids.data());

    // paint encoding for the point number, 
    float point_rate = (float)point_count_ / max_point_count_;
    int seg_num = (SEG_PER_CIRCLE - 1) * point_rate;
    std::vector<vtkIdType > num_inner_ids, num_outer_ids;
    for (int j = 0; j <= seg_num; ++j) {
        float end_arc = -1 * j * 3.14159 * 2 / SEG_PER_CIRCLE + 3.14159 * 0.5;
        float cos_value = cos(end_arc);
        float sin_value = sin(end_arc);
        float x = node_radius_ * cos_value * 1.05;
		float y = node_radius_ * sin_value * 1.05;

        vtkIdType id_one = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.0000001);
		glyph_colors->InsertNextTuple4(200, 200, 200, 255);
        num_inner_ids.push_back(id_one);

        x = node_radius_ * cos_value * 1.18;
		y = node_radius_ * sin_value * 1.18;
        vtkIdType id_two = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.0000001);
		glyph_colors->InsertNextTuple4(200, 200, 200, 255);
        num_outer_ids.push_back(id_two);
    }

    vtkIdType cell_ids[3];
    for (int j = 0; j < num_outer_ids.size() - 1; ++j) {
        cell_ids[0] = num_outer_ids[j];
        cell_ids[1] = num_outer_ids[j + 1];
        cell_ids[2] = num_inner_ids[j];
        glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

        cell_ids[0] = num_outer_ids[j + 1];
        cell_ids[1] = num_inner_ids[j + 1];
        cell_ids[2] = num_inner_ids[j];
        glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
    }

    vector<vtkIdType> center_cirlce_ids;
	for (int j = 0; j <= SEG_PER_PIE; ++j) {
		float end_arc = j * 3.14159 * 2 / SEG_PER_PIE;
		float x = node_radius_ * 0.05 * cos(end_arc);
		float y = node_radius_ * 0.05 * sin(end_arc);

		center_cirlce_ids.push_back(glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.001));
		glyph_colors->InsertNextTuple4(0, 0, 0, 255);
	}
	glyph_poly_->InsertNextCell(VTK_POLYGON, SEG_PER_PIE + 1, center_cirlce_ids.data());

    // paint glyph
    std::vector<vtkIdType> outer_band_ids;
    std::vector<vtkIdType> inner_band_ids;
    std::vector<vtkIdType> gray_region_ids;

    outer_band_ids.resize(2 * SEG_PER_PIE);
    inner_band_ids.resize(2 * SEG_PER_PIE);
    gray_region_ids.resize(SEG_PER_PIE + 1);

    for (int i = 0; i < var_num; ++i) {
        float begin_arc = i * PIE_VAL * 2 / var_num;
        float end_arc = (i + 1) * PIE_VAL * 2 / var_num;
        float step_arc = (end_arc - begin_arc) / (SEG_PER_PIE - 1);

        QColor pie_color = colors_[i];

        float temp_arc = begin_arc;
        for (int j = 0; j < SEG_PER_PIE; ++j) {
            float cos_value = cos(temp_arc);
            float sin_value = sin(temp_arc);
            float std_dev = std_devs_[i] < 0.03 ? 0.03 : std_devs_[i];

            // insert the outer polygon
            float temp_radius = node_radius_ * (means_[i] + std_dev) * 0.9 + node_radius_ * 0.1;
            if (temp_radius > node_radius_) temp_radius = node_radius_;
            float x = temp_radius * cos_value;
            float y = temp_radius * sin_value;
            vtkIdType id_one = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0);
            glyph_colors->InsertNextTuple4(pie_color.red(), pie_color.green(), pie_color.blue(), 64);
            outer_band_ids[2 * j + 1] = id_one;

            // insert the inner polygon
            temp_radius = node_radius_ * (means_[i] - std_dev) * 0.9 + node_radius_ * 0.1;
            if (temp_radius < 0.3 * node_radius_) temp_radius = 0.3 * node_radius_;
            x = temp_radius * cos_value;
            y = temp_radius * sin_value;

            vtkIdType id_two = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0);
            glyph_colors->InsertNextTuple4(pie_color.red(), pie_color.green(), pie_color.blue(), 64);
            inner_band_ids[2 * j + 1] = id_two;

            // insert the gray polygon
            vtkIdType id_four = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0);
            glyph_colors->InsertNextTuple4(240, 240, 240, 255);
            gray_region_ids[j] = id_four;

            // insert the average value
            temp_radius = node_radius_ * means_[i] * 0.9 + node_radius_ * 0.1;
            x = temp_radius * cos_value;
            y = temp_radius * sin_value;
            vtkIdType id_three = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0);
            glyph_colors->InsertNextTuple4(pie_color.red(), pie_color.green(), pie_color.blue(), 255);
            outer_band_ids[2 * j] = id_three;
            inner_band_ids[2 * j] = id_three;

            temp_arc += step_arc;
        }

        vtkIdType center_id = glyph_points->InsertNextPoint(center_x_, center_y_, 0);
        glyph_colors->InsertNextTuple4(240, 240, 240, 255);
        gray_region_ids[SEG_PER_PIE] = center_id;

        glyph_poly_->InsertNextCell(VTK_POLYGON, (int)gray_region_ids.size(), gray_region_ids.data());
        glyph_poly_->InsertNextCell(VTK_TRIANGLE_STRIP, (int)outer_band_ids.size(), outer_band_ids.data());
        glyph_poly_->InsertNextCell(VTK_TRIANGLE_STRIP, (int)inner_band_ids.size(), inner_band_ids.data());
    }

    std::vector<vtkIdType> range_outer_ids, range_inner_ids, saliency_ids;
    for (int j = 0; j < SEG_PER_CIRCLE; ++j) {
        float end_arc = -1 * j * 3.14159 * 2 / (SEG_PER_CIRCLE - 1);
        float cos_value = cos(end_arc);
        float sin_value = sin(end_arc);

        float x = node_radius_ * cos_value;
        float y = node_radius_ * sin_value;
        vtkIdType id_three = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0.0001);
		glyph_colors->InsertNextTuple4(150, 150, 150, 255);
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
        float begin_arc = highlight_index * PIE_VAL * 2 / var_num;
        float end_arc = (highlight_index + 1) * PIE_VAL * 2 / var_num;
        float step_arc = (end_arc - begin_arc) / (SEG_PER_PIE - 1);

		vector<int> comp_ids;
		comp_ids.resize(SEG_PER_PIE);

		float temp_arc = begin_arc;
		for (int k = 0; k < SEG_PER_PIE; ++k) {
			float cos_value = cos(temp_arc);
			float sin_value = sin(temp_arc);

			// insert the average value
			float temp_radius = node_radius_ * highlight_val * 0.9 + node_radius_ * 0.1;
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
		    comp_ids[k] = highlight_points->InsertNextPoint(center_x_ + x, center_y_ + y, 1e-3);
            highlight_colors->InsertNextTuple4(255, 0, 0, 255);

		    temp_arc += step_arc;
	    }

        highlight_poly_->InsertNextCell(VTK_POLY_LINE, comp_ids.size(), comp_ids.data());
    }
}

void RadarBandWidgetNoSaliency::BuildSmallGlyph() {
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

	vtkIdType center_id = glyph_points->InsertNextPoint(center_x_, center_y_, 0);
	vtkIdType pre_id = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0);

    vtkIdType cell_ids[3];
	for (int k = 1; k <= 20; ++k) {
		end_arc = (float)k / 20 * 3.14159 * 2;
		x = temp_radius * cos(end_arc);
		y = temp_radius * sin(end_arc);
		vtkIdType current_id = glyph_points->InsertNextPoint(center_x_ + x, center_y_ + y, 0);

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

void RadarBandWidgetNoSaliency::OnMouseMove() {
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

	    int highlight_index = (int)(degree / (2 * 3.1416) * glyph_object_->var_num());
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

void RadarBandWidgetNoSaliency::OnLeftButtonDown() {
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

void RadarBandWidgetNoSaliency::OnLeftButtonUp() {
    this->glyph_state_ = GlyphWidget::NORMAL;
}

void RadarBandWidgetNoSaliency::OnRightButtonDown() {

}

void RadarBandWidgetNoSaliency::OnRightButtonUp() {

}