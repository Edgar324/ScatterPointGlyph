/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef MEAN_STD_WIDGET_H_
#define MEAN_STD_WIDGET_H_

#include "glyph_widget.h"

class MeanStdWidget : public GlyphWidget
{
public:
    static MeanStdWidget* New();

    vtkTypeMacro(MeanStdWidget, GlyphWidget);

    virtual void SetEnabled(int enabled);
    void SetRenderingMode(int mode);

protected:
    MeanStdWidget();
    ~MeanStdWidget();

    virtual void BuildRepresentation();
    virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();

private:
    int rendering_mode_ = 0;

    vtkActor* glyph_actor_ = NULL;
    vtkPolyDataMapper* glyph_data_mapper_ = NULL;
    vtkPolyData* glyph_poly_ = NULL;

    vtkActor* highlight_actor_ = NULL;
    vtkPolyDataMapper* highlight_data_mapper_ = NULL;
    vtkPolyData* highlight_poly_ = NULL;

    void BuildLargeVarGlyph();
    void BuildLargeMeanGlyph();
    void BuildSmallGlyph();
    void MakeCubic(vector<QPointF>& points, int point_num, vector<QPointF>& splines);
};

#endif