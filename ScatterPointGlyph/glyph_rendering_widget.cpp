/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "glyph_rendering_widget.h"
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <QVBoxLayout>
#include "glyph_object.h"
#include "glyph_dataset.h"
#include "radar_band_widget.h"
#include "regular_radar_widget.h"
#include "radar_band_widget_no_saliency.h"
#include "qvtk_rendering_widget.h"
#include "point_rendering_widget.h"
#include "arrow_widget.h"
#include "density_rendering_widget.h"
#include "map_rendering_widget.h"
#include "indicator_widget.h"


GlyphRenderingWidget::GlyphRenderingWidget() {
    rendering_widget_ = new QVtkRenderingWidget;
	rendering_widget_->setFocusPolicy(Qt::StrongFocus);

    connect(rendering_widget_, SIGNAL(ViewUpdated()), this, SIGNAL(ViewportUpdated()));
    connect(rendering_widget_, SIGNAL(SelectionChanged()), this, SIGNAL(BrushSelectionChanged()));

    map_widget_ = MapRenderingWidget::New();
    map_widget_->SetInteractor(rendering_widget_->GetInteractor());

    point_rendering_widget_ = PointRenderingWidget::New();
    point_rendering_widget_->SetInteractor(rendering_widget_->GetInteractor());

    arrow_widget_ = ArrowWidget::New();
    arrow_widget_->SetInteractor(rendering_widget_->GetInteractor());

    density_widget_ = DensityRenderingWidget::New();
    density_widget_->SetInteractor(rendering_widget_->GetInteractor());

    indicator_widget_ = IndicatorWidget::New();
    indicator_widget_->SetInteractor(rendering_widget_->GetInteractor());

    QVBoxLayout* main_layout = new QVBoxLayout;
    main_layout->addWidget(rendering_widget_);
    this->setLayout(main_layout);
}

GlyphRenderingWidget::~GlyphRenderingWidget() {

}

void GlyphRenderingWidget::SetMvDataset(MultivariateDataset* dataset) {
    this->mv_dataset_ = dataset;

    // initialize the point rendering layer
    point_rendering_widget_->SetEnabled(false);
    point_rendering_widget_->SetData(mv_dataset_);
    point_rendering_widget_->SetEnabled(true);

    arrow_widget_->SetEnabled(true);

    indicator_widget_->SetRenderer(rendering_widget_->indicator_renderer());
    indicator_widget_->SetData(0);
    indicator_widget_->SetEnabled(true);

    rendering_widget_->main_renderer()->ResetCamera();
    rendering_widget_->update();
}

void GlyphRenderingWidget::SetData(GlyphDataset* glyph_dataset) {
    this->glyph_dataset_ = glyph_dataset;
    selected_objects_.clear();

    connect(this->glyph_dataset_, SIGNAL(DataUpdated()), this, SLOT(OnDataUpdated()));

    this->UpdateView();
}

void GlyphRenderingWidget::SetDensityMapVisibility(bool visible) {
    this->density_widget_->SetEnabled(visible);
    rendering_widget_->update();
}

void GlyphRenderingWidget::SetPointMapVisibility(bool visible) {
    this->point_rendering_widget_->SetEnabled(visible);
    rendering_widget_->update();
}

void GlyphRenderingWidget::SetGlyphVisibility(bool visible) {
    is_glyph_visible_ = visible;

    vector<GlyphObject*> objects;
    this->glyph_dataset_->GetAllGlyphObjects(objects);
    for (int i = 0; i < objects.size(); ++i) {
        this->glyph_widgets_[i]->SetEnabled(visible);
    }
    rendering_widget_->update();
}

void GlyphRenderingWidget::SetGeoMapVisibility(bool visible) {
    map_widget_->SetEnabled(visible);
    rendering_widget_->update();
}

void GlyphRenderingWidget::SetWidgetState(WidgetState state) {
    this->widget_state_ = state;

    // Clear all widget state
    list<GlyphObject*>::iterator iter = selected_objects_.begin();
    while (iter != selected_objects_.end()) {
        (*iter)->set_is_selected(false);
        iter++;
    }
    selected_objects_.clear();

    if (this->widget_state_ == BRUSH_SELECTION) {
        this->rendering_widget_->SetBrushEnabled(true);
    } else {
        this->rendering_widget_->SetBrushEnabled(false);
    }
}

void GlyphRenderingWidget::SetSelectedPointIds(vector<vector<int>>& ids) {
    this->point_rendering_widget_->SetSelectedIds(ids);
}

void GlyphRenderingWidget::SetColorMapping(QString name, vector<float>& values, float min_val, float max_val) {
    this->point_rendering_widget_->SetColorMappingOn(name, values, min_val, max_val);
}

void GlyphRenderingWidget::SetColorMappingOff() {
    this->point_rendering_widget_->SetColorMappingOff();
}

void GlyphRenderingWidget::SetPointDensityInfo(vector<vector<float>>& pos, vector<int>& cluster_index) {
    this->density_widget_->SetData(pos, cluster_index);
}

void GlyphRenderingWidget::GetViewPort(float& left, float& right, float& bottom, float& top) {
    this->rendering_widget_->GetViewport(left, right, bottom, top);
}

void GlyphRenderingWidget::GetSelectedClusters(vector<int>& cluster_ids) {
    cluster_ids.clear();
    list<GlyphObject*>::iterator iter = selected_objects_.begin();
    while (iter != selected_objects_.end()) {
        cluster_ids.push_back((*iter)->cluster_id());
        iter++;
    }
}

void GlyphRenderingWidget::GetBrushingPath(vector<float>& path) {
    this->rendering_widget_->GetBrushingPath(path);
}

void GlyphRenderingWidget::UpdateView() {
    this->BuildRepresentation();
    this->UpdateArrowWidget();
    indicator_widget_->SetData(this->glyph_dataset_->max_point_count());
}

void GlyphRenderingWidget::BuildRepresentation() {
    // Update rendering widget
    vector<GlyphObject*> objects;
    this->glyph_dataset_->GetAllGlyphObjects(objects);

    if (this->glyph_widgets_.size() < objects.size()) {
        int cur_size = this->glyph_widgets_.size();
        for (int i = cur_size; i < objects.size(); ++i) {
            GlyphWidget* widget = RadarBandWidget::New();
            this->glyph_widgets_.push_back(widget);
            widget->SetInteractor(rendering_widget_->GetInteractor());
            widget->SetRenderingWidget(this);
        }
    }

    for (int i = 0; i < objects.size(); ++i) {
        this->glyph_widgets_[i]->SetData(objects[i]);
        this->glyph_widgets_[i]->SetEnabled(is_glyph_visible_);
    }

    for (int i = objects.size(); i < this->glyph_widgets_.size(); ++i) {
        this->glyph_widgets_[i]->SetData(NULL);
        this->glyph_widgets_[i]->SetEnabled(false);
    }

    this->rendering_widget_->update();

    // update point rendering view
    /*float left, right, bottom, top;
    this->rendering_widget_->GetViewport(left, right, bottom, top);
    this->point_rendering_widget_->SetViewpot(left, right, bottom, top);*/
}

void GlyphRenderingWidget::ClearRepresentation() {
    for (int i = 0; i < glyph_widgets_.size(); ++i) {
        glyph_widgets_[i]->SetEnabled(false);
        glyph_widgets_[i]->Delete();
    }
    glyph_widgets_.clear();
}

float GlyphRenderingWidget::GetDistancePerPixel() {
    return rendering_widget_->GetDistancePerPixel();
}

void GlyphRenderingWidget::OnDataUpdated() {
    this->selected_objects_.clear();
    this->UpdateView();
}

void GlyphRenderingWidget::HighlightModified(GlyphObject* object) {
    int temp_index = object->highlight_index();
    float mean_val = object->highlight_val();
    float std_val = object->std_devs()[temp_index];

    if (temp_index != -1) {
        QString tip;
        tip += object->names()[temp_index] + QString(" Mean: ") + QString::number(mean_val) + QString("\t\n");
        tip += object->names()[temp_index] + QString(" SD: ") + QString::number(std_val) + QString("\t\n");
        tip += QString("Point Count: ") + QString::number(object->point_count());
        this->rendering_widget_->ShowTooltip(tip);
    } else {
        this->rendering_widget_->HideTooltip();
    }

    if (temp_index != highlight_index_) {
        highlight_index_ = temp_index;
        highlight_val_ = mean_val;

        vector<GlyphObject*> glyph_objects;
        glyph_dataset_->GetAllGlyphObjects(glyph_objects);
        for (int i = 0; i < glyph_objects.size(); ++i) {
            glyph_objects[i]->set_highlight(highlight_index_, highlight_val_);
        }
    }
}

void GlyphRenderingWidget::SelectionModified(GlyphObject* object) {
    list<GlyphObject*>::iterator iter = selected_objects_.begin();
    while (iter != selected_objects_.end() 
        && (*iter)->cluster_id() != object->cluster_id()) 
        iter++;

    if (object->is_selected()) {
        if (iter == selected_objects_.end()) selected_objects_.push_back(object);
    } else {
        if (iter != selected_objects_.end()) selected_objects_.erase(iter);
    }

    this->UpdateArrowWidget();

    emit SelectionChanged();
}

void GlyphRenderingWidget::UpdateArrowWidget() {
    vector<float> point_pos;
    float glyph_radius = 0;
    list<GlyphObject*>::iterator iter = selected_objects_.begin();
    if (iter != selected_objects_.end()) {
        float pre_x = (*iter)->center_x();
        float pre_y = (*iter)->center_y();
        iter++;
        while (iter != selected_objects_.end()) {
            point_pos.push_back(pre_x);
            point_pos.push_back(pre_y);
            point_pos.push_back((*iter)->center_x());
            point_pos.push_back((*iter)->center_y());
            pre_x = (*iter)->center_x();
            pre_y = (*iter)->center_y();
            glyph_radius = (*iter)->node_radius();
            iter++;
        }
    }
    arrow_widget_->SetData(point_pos, glyph_radius);
}