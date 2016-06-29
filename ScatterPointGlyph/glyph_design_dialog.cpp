/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "glyph_design_dialog.h"
#include <QVTKWidget.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkTextActor3D.h>
#include "vtkTextProperty.h"

GlyphDesignDialog::GlyphDesignDialog()
{
    ui_.setupUi(this);

    vtkwidget_ = new QVTKWidget;
    main_renderer_ = vtkRenderer::New();
	vtkwidget_->GetRenderWindow()->AddRenderer(main_renderer_);
	main_renderer_->SetViewport(0.0, 0.0, 1.0, 1.0);
	main_renderer_->SetBackground(1.0, 1.0, 1.0);

    glyph_actor_ = vtkActor::New();
    glyph_mapper_ = vtkPolyDataMapper::New();
    glyph_poly_ = vtkPolyData::New();
    glyph_mapper_->SetInputData(glyph_poly_);
    glyph_actor_->SetMapper(glyph_mapper_);
    main_renderer_->AddActor(glyph_actor_);

    QHBoxLayout* glyph_layout = new QHBoxLayout;
    glyph_layout->addWidget(vtkwidget_);
    ui_.glyph_groupbox->setLayout(glyph_layout);
}

GlyphDesignDialog::~GlyphDesignDialog()
{
    for (int i = 0; i < variable_checkbox_.size(); ++i) {
        delete variable_checkbox_[i];
    }
    variable_checkbox_.clear();
}

void GlyphDesignDialog::SetVariables(std::vector<QString >& names)
{
    mean_values_.clear();
    msd_values_.clear();
    QGridLayout* gridlayout = new QGridLayout;
    for (int i = 0; i < names.size(); ++i) {
        QCheckBox* check_box = new QCheckBox(names[i]);
        check_box->setChecked(true);
        variable_checkbox_.push_back(check_box);
        gridlayout->addWidget(check_box, i / 2, i % 2, 1, 1);

        connect(check_box, SIGNAL(stateChanged(int)), this, SLOT(OnVariableSelectionChanged()));

        mean_values_.push_back((float)rand() / RAND_MAX);
        float msd = (float)rand() / RAND_MAX;
        while (msd > mean_values_[i]) msd /= 5;
        msd_values_.push_back(msd);
    }
    ui_.variable_groupbox->setLayout(gridlayout);

    UpdateGlyph();
}

void GlyphDesignDialog::GetSelectedVariables(std::vector<int>& index)
{
    index.clear();
    for (int i = 0; i < variable_checkbox_.size(); ++i)
        if (variable_checkbox_[i]->isChecked()) index.push_back(i);
}

void GlyphDesignDialog::UpdateGlyph() {
    glyph_poly_->Initialize();

    std::vector<int> selection_index;
    this->GetSelectedVariables(selection_index);

    float node_center_x = 0.0;
    float node_center_y = 0.0;
    float node_radius = 1.0;
    float saliency = 0.5;
    int var_num = selection_index.size();

    vtkPoints* points = vtkPoints::New();
	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
	colors->SetNumberOfComponents(4);

	glyph_poly_->SetPoints(points);
	glyph_poly_->GetPointData()->SetScalars(colors);

	vtkCellArray* poly_array = vtkCellArray::New();
	glyph_poly_->SetPolys(poly_array);

	vtkCellArray* line_array = vtkCellArray::New();
	glyph_poly_->SetLines(line_array);

	vtkCellArray* strip_array = vtkCellArray::New();
	glyph_poly_->SetStrips(strip_array);

    int seg_per_circle = 50;
	// insert background
	std::vector<vtkIdType > background_ids;
	for (int j = 0; j <= seg_per_circle; ++j) {
		float end_arc = j * 3.14159 * 2 / seg_per_circle;
		float x = node_radius * cos(end_arc);
		float y = node_radius * sin(end_arc);

		background_ids.push_back(points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.001));
		colors->InsertNextTuple4(255, 255, 255, 10);
	}
	glyph_poly_->InsertNextCell(VTK_POLYGON, seg_per_circle + 1, background_ids.data());

	std::vector<vtkIdType > center_cirlce_ids;
	for (int j = 0; j <= seg_per_circle; ++j) {
		float end_arc = j * 3.14159 * 2 / seg_per_circle;
		float x = node_radius * 0.05 * cos(end_arc);
		float y = node_radius * 0.05 * sin(end_arc);

		center_cirlce_ids.push_back(points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.003));
		colors->InsertNextTuple4(0, 0, 0, 255);
	}
	glyph_poly_->InsertNextCell(VTK_POLYGON, seg_per_circle + 1, center_cirlce_ids.data());

	std::vector<vtkIdType > radius_circle_ids;
	for (int j = 0; j <= seg_per_circle; ++j) {
		float end_arc = j * 3.14159 * 2 / seg_per_circle;
		float x = node_radius * cos(end_arc);
		float y = node_radius * sin(end_arc);

		radius_circle_ids.push_back(points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.001));
        float gray_value = 250 * (1.0 - exp(-3 * (1.0 - saliency)));
		colors->InsertNextTuple4(gray_value, gray_value, gray_value, 255);
	}

	vtkIdType circle_ids[2];
	for (int j = 0; j < seg_per_circle; ++j) {
		circle_ids[0] = radius_circle_ids[j];
		circle_ids[1] = radius_circle_ids[j + 1];
		glyph_poly_->InsertNextCell(VTK_LINE, 2, circle_ids);
	}

    // paint encoding for the point number
    float point_rate = 0.5;
    int seg_num = seg_per_circle * point_rate;
    std::vector<vtkIdType > inner_ids, outer_ids;
    for (int j = 0; j <= seg_num; ++j) {
        float end_arc = -1 * j * 3.14159 * 2 / seg_per_circle + 3.14159 * 0.5;
        float x = node_radius * cos(end_arc) * 1.05;
		float y = node_radius * sin(end_arc) * 1.05;

        vtkIdType id_one = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
		colors->InsertNextTuple4(128, 128, 128, 128);
        inner_ids.push_back(id_one);

        x = node_radius * cos(end_arc) * 1.2;
		y = node_radius * sin(end_arc) * 1.2;

        vtkIdType id_two = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
		colors->InsertNextTuple4(128, 128, 128, 128);
        outer_ids.push_back(id_two);
    }

    vtkIdType cell_ids[5];
    for (int j = 0; j < outer_ids.size() - 1; ++j) {
        cell_ids[0] = outer_ids[j];
        cell_ids[1] = outer_ids[j + 1];
        cell_ids[2] = inner_ids[j];
        glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

        cell_ids[0] = outer_ids[j + 1];
        cell_ids[1] = inner_ids[j + 1];
        cell_ids[2] = inner_ids[j];
        glyph_poly_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
    }

	// paint variance
	std::vector<vtkIdType > var_point_ids1;
	std::vector<vtkIdType > var_point_ids2;
	std::vector<vtkIdType > gray_region_ids;
	int seg_per_pie = 30;
	var_point_ids1.resize(2 * seg_per_pie);
	var_point_ids2.resize(2 * seg_per_pie);
	gray_region_ids.resize(seg_per_pie + 1);

	vtkIdType line_ids[2];
	float alpha = 1.0;

	for (int j = 0; j < var_num; ++j) {
		float begin_arc = j * 3.14159 * 2 / var_num;
		float end_arc = (j + 1) * 3.14159 * 2 / var_num;
		float step_arc = (end_arc - begin_arc) / (seg_per_pie - 1);

		int var_index = selection_index[j];

		float temp_arc = begin_arc;
		for (int k = 0; k < seg_per_pie; ++k) {
			float cos_value = cos(temp_arc);
			float sin_value = sin(temp_arc);

			// insert the outer polygon
			float temp_radius = node_radius * (mean_values_[var_index] + msd_values_[var_index]) * 0.8 + node_radius * 0.1;
			if (temp_radius > node_radius) temp_radius = node_radius;
			float x = temp_radius * cos_value;
			float y = temp_radius * sin_value;
			vtkIdType id_one = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
			colors->InsertNextTuple4(70, 255, 70, 128 * alpha);
			var_point_ids1[2 * k + 1] = id_one;

			// insert the inner polygon
			temp_radius = node_radius * (mean_values_[var_index] - msd_values_[var_index]) * 0.8 + node_radius * 0.1;
			if (temp_radius < 0) temp_radius = 0;
			x = temp_radius * cos_value;
			y = temp_radius * sin_value;

			vtkIdType id_two = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
			colors->InsertNextTuple4(70, 255, 70, 128 * alpha);
			var_point_ids2[2 * k + 1] = id_two;

			// insert the gray polygon
			vtkIdType id_four = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
			colors->InsertNextTuple4(230, 230, 230, 255 * alpha);
			gray_region_ids[k] = id_four;

			// insert the average value
			temp_radius = node_radius * mean_values_[var_index] * 0.8 + node_radius * 0.1;
			x = temp_radius * cos_value;
			y = temp_radius * sin_value;
			vtkIdType id_three = points->InsertNextPoint(node_center_x + x, node_center_y + y, 0.002);
			colors->InsertNextTuple4(70, 255, 70, 255 * alpha);
			var_point_ids1[2 * k] = id_three;
			var_point_ids2[2 * k] = id_three;

			temp_arc += step_arc;
		}
				
		vtkIdType center_id = points->InsertNextPoint(node_center_x, node_center_y, 0.002);
		colors->InsertNextTuple4(230, 230, 230, 255 * alpha);
		gray_region_ids[seg_per_pie] = center_id;

		glyph_poly_->InsertNextCell(VTK_POLYGON, (int)gray_region_ids.size(), gray_region_ids.data());
		glyph_poly_->InsertNextCell(VTK_TRIANGLE_STRIP, (int)var_point_ids1.size(), var_point_ids1.data());
		glyph_poly_->InsertNextCell(VTK_TRIANGLE_STRIP, (int)var_point_ids2.size(), var_point_ids2.data());

		for (int k = 0; k < seg_per_pie - 1; ++k) {
			line_ids[0] = var_point_ids1[2 * k];
			line_ids[1] = var_point_ids1[2 * (k + 1)];
			glyph_poly_->InsertNextCell(VTK_LINE, 2, line_ids);
		}
	}

    // update text actors
    while (text_actors_.size() > 0) {
        main_renderer_->RemoveActor(text_actors_[text_actors_.size() - 1]);
        text_actors_[text_actors_.size() - 1]->Delete();
        text_actors_.resize(text_actors_.size() - 1);
    }
    for (int j = 0; j < selection_index.size(); ++j) {
        vtkTextActor3D* actor = vtkTextActor3D::New();
		actor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
		actor->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
		actor->GetTextProperty()->SetBold(true);
		actor->GetTextProperty()->SetFontFamilyToArial();
		main_renderer_->AddActor(actor);
		text_actors_.push_back(actor);
    }

    float accu_right = 0, accu_left = 0;
    float text_height = 2 * node_radius / (var_num / 2 + 1);
    for (int i = 0; i < var_num; ++i) {
  //      // insert indicator line
  //      float center_arc = (i + 0.5) * 3.14159 * 2 / var_num;
  //      float text_pos[2], begin_pos[2], arc_pos[2], end_pos[2];
  //      if ((i + 0.5) / var_num > 0.25 && ((i + 0.5) / var_num < 0.75)) {
  //          text_pos[0] = node_center_x - 3 * node_radius;
  //          text_pos[1] = node_center_y + node_radius - accu_left * text_height;
  //          
  //          begin_pos[0] = node_center_x + node_radius * cos(center_arc) * 0.5;
  //          begin_pos[1] = node_center_y + node_radius * sin(center_arc) * 0.5;

  //          arc_pos[0] = node_center_x - node_radius * cos(center_arc) * 1.5;
  //          arc_pos[1] = node_center_y + node_radius * sin(center_arc) * 1.5;

  //          end_pos[0] = node_center_x - node_radius * 2;
  //          end_pos[1] = arc_pos[1];

  //          accu_left++;
  //      }
  //      else {
  //          text_pos[0] = node_center_x + 2 * node_radius;
  //          text_pos[1] = node_center_y + node_radius - accu_right * text_height;

  //          begin_pos[0] = node_center_x + node_radius * cos(center_arc) * 0.5;
  //          begin_pos[1] = node_center_y + node_radius * sin(center_arc) * 0.5;

  //          arc_pos[0] = node_center_x - node_radius * cos(center_arc) * 1.5;
  //          arc_pos[1] = node_center_y + node_radius * sin(center_arc) * 1.5;

  //          end_pos[0] = node_center_x + node_radius * 2;
  //          end_pos[1] = arc_pos[1];

  //          accu_right++;
  //      }
  //      vtkIdType line_ids[2];
  //      line_ids[0] = points->InsertNextPoint(begin_pos[0], begin_pos[1], 0.002);
		//colors->InsertNextTuple4(80, 80, 80, 255 * alpha);

  //      line_ids[1] = points->InsertNextPoint(arc_pos[0], arc_pos[1], 0.002);
		//colors->InsertNextTuple4(80, 80, 80, 255 * alpha);

  //      glyph_poly_->InsertNextCell(VTK_LINE, 2, line_ids);

  //      line_ids[0] = points->InsertNextPoint(end_pos[0], end_pos[1], 0.002);
		//colors->InsertNextTuple4(80, 80, 80, 255 * alpha);

  //      glyph_poly_->InsertNextCell(VTK_LINE, 2, line_ids);

        float text_pos[2];
        float center_arc = (i + 0.5) * 3.14159 * 2 / var_num;
        text_pos[0] = node_center_x + node_radius * cos(center_arc) * 1.5;
        text_pos[1] = node_center_y + node_radius * sin(center_arc) * 1.5;

		text_actors_[i]->SetInput(variable_checkbox_[selection_index[i]]->text().toLocal8Bit().data());
		text_actors_[i]->SetPosition(text_pos[0], text_pos[1], 0.003);
		text_actors_[i]->GetTextProperty()->SetFontSize(20);
		float scale = node_radius / 180;
		text_actors_[i]->SetScale(scale, scale, scale);
        text_actors_[i]->RotateWXYZ((i + 0.5) / var_num * 360, 0, 0, 1);
		text_actors_[i]->Modified();
    }

    glyph_actor_->Modified();
    vtkwidget_->update();
}

void GlyphDesignDialog::OnVariableSelectionChanged()
{
    this->UpdateGlyph();
}
