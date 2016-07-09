/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "glyph_object.h"
#include "glyph_widget.h"

GlyphObject::GlyphObject(int cluster_id, vector<QString>& names, vector<QColor>& colors,
    vector<float>& means, vector<float>& std_devs,
    float saliency, int point_count, int max_point_count,
    float node_radius, float center_x, float center_y, bool is_expandable)
    : cluster_id_(cluster_id), names_(names), colors_(colors), means_(means), std_devs_(std_devs),
    saliency_(saliency), point_count_(point_count), max_point_count_(max_point_count),
    node_radius_(node_radius), center_x_(center_x), center_y_(center_y), is_expandable_(is_expandable),
    highlight_index_(-1), is_selected_(false) {

}

GlyphObject::~GlyphObject() {

}

void GlyphObject::set_highlight(int index, float val) {
    if (index != highlight_index_ || val != highlight_val_) {
        highlight_index_ = index;
        highlight_val_ = val;
        this->glyph_widget_->Modified();
    }
}

void GlyphObject::set_is_selected(bool selected) {
    if (is_selected_ != selected) {
        is_selected_ = selected;
        this->glyph_widget_->Modified();
    }
}

void GlyphObject::set_glyph_widget(GlyphWidget* widget) {
    this->glyph_widget_ = widget;
}
