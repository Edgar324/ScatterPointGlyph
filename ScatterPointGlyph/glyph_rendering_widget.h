/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef GLYPH_RENDERING_WIDGET_H_
#define GLYPH_RENDERING_WIDGET_H_

#include <vector>
#include <list>
using namespace std;
#include <QWidget>
#include "glyph_dataset.h"

class QVtkRenderingWidget;
class ScatterPointDataset;
class GlyphWidget;
class ArrowWidget;
class PointRenderingWidget;
class DensityRenderingWidget;

class GlyphRenderingWidget : public QWidget
{
    Q_OBJECT

public:
    GlyphRenderingWidget();
    ~GlyphRenderingWidget();

    enum WidgetState {
        NORMAL = 0x0,
        MOUSE_SELECTION,
        BRUSH_SELECTION
    };

    void SetDensityMapVisibility(bool visible);
    void SetPointMapVisibility(bool visible);
    void SetGlyphVisibility(bool visible);
    void SetGeoMapVisibility(bool visible);

    void SetPointData(ScatterPointDataset* point_dataset);
    void SetData(GlyphDataset* glyph_dataset);
    // for point rendering widget
    void SetSelectedPointIds(vector<vector<int>>& ids);
    void SetColorMapping(QString name, vector<float>& values, float min_val, float max_val);
    void SetColorMappingOff();
    // for density rendering widget
    void SetPointDensityInfo(vector<vector<float>>& pos, vector<int>& cluster_index);

    void SetWidgetState(WidgetState state);

    void GetViewPort(float& left, float& right, float& bottom, float& top);
    float GetDistancePerPixel();

    void GetSelectedClusters(vector<int>& cluster_ids);
    void GetBrushingPath(vector<float>& path);

    void HighlightModified(GlyphObject* object);
    void SelectionModified(GlyphObject* object);

signals:
    void ViewportUpdated();
    void SelectionChanged();
    void BrushSelectionChanged();

protected:

private:
    ScatterPointDataset* point_dataset_ = NULL;
    GlyphDataset* glyph_dataset_ = NULL;

    bool is_glyph_visible_ = true;
    QVtkRenderingWidget* rendering_widget_ = NULL;
    vector<GlyphWidget*> glyph_widgets_;

    PointRenderingWidget* point_rendering_widget_ = NULL;

    ArrowWidget* arrow_widget_;

    DensityRenderingWidget* density_widget_;

    WidgetState widget_state_ = NORMAL;
    bool is_multiple_selection = false;

    list<GlyphObject*> selected_objects_;

    int highlight_index_ = -1;
    float highlight_val_ = -1;

    void UpdateView();
    void UpdateArrowWidget();
    void BuildRepresentation();
    void ClearRepresentation();

private slots:
    void OnDataUpdated();
};

#endif