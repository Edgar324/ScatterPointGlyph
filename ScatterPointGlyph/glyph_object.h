/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef GLYPH_OBJECT_H_
#define GLYPH_OBJECT_H_

#include <vector>
#include <map>
using namespace std;
#include <QString>
#include <QColor>

class GlyphWidget;

class GlyphObject : public QObject
{
    Q_OBJECT

public:
    GlyphObject(int cluster_id, vector<QString>& names, vector<QColor>& colors,
        vector<float>& means, vector<float>& std_devs, 
        float saliency, int point_count, int max_point_count,
        float node_radius, float center_x, float center_y, bool is_expandable);
    ~GlyphObject();

    float node_radius() { return node_radius_; }
    float center_x() { return center_x_; }
    float center_y() { return center_y_; }

    int cluster_id() { return cluster_id_; }
    int point_count() { return point_count_; }
    int max_point_count() { return max_point_count_; }
    float saliency() { return saliency_; }
    int var_num() { return means_.size(); }
    vector<float>& means() { return means_; }
    vector<float>& std_devs() { return std_devs_; }
    vector<QString>& names() { return names_; }
    vector<QColor>& colors() { return colors_; }
    int highlight_index() { return highlight_index_; }
    float highlight_val() { return highlight_val_; }
    bool is_selected() { return is_selected_; }
    bool is_expandable() { return is_expandable_; }

    void set_highlight(int index, float val);
    void set_is_selected(bool selected);

    void set_glyph_widget(GlyphWidget* widget);

signals:
    void HighlightChanged();
    void SelectionStatusChanged();

private:
    // glyph representation
    float node_radius_;
    float center_x_, center_y_;

    // data
    int cluster_id_;
    int point_count_; 
    int max_point_count_;
    float saliency_;
    vector<QString> names_;
    vector<QColor> colors_;
    vector<float> means_;
    vector<float> std_devs_;

    int highlight_index_ = 0;
    float highlight_val_ = 0;
    bool is_selected_ = false;
    bool is_expandable_ = false;

    GlyphWidget* glyph_widget_ = NULL;
};

#endif