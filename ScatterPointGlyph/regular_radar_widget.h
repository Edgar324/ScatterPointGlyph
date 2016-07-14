/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef REGULAR_RADAR_WIDGET_H_
#define REGULAR_RADAR_WIDGET_H_

#include "glyph_widget.h"

class RegularRadarWidget : public GlyphWidget
{
public:
    static RegularRadarWidget* New();

    vtkTypeMacro(RegularRadarWidget, GlyphWidget);

    virtual void SetEnabled(int enabled);

protected:
    RegularRadarWidget();
    ~RegularRadarWidget();

    virtual void BuildRepresentation();
    virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();

private:
    vtkActor* glyph_actor_ = NULL;
    vtkPolyDataMapper* glyph_data_mapper_ = NULL;
    vtkPolyData* glyph_poly_ = NULL;

    vtkActor* highlight_actor_ = NULL;
    vtkPolyDataMapper* highlight_data_mapper_ = NULL;
    vtkPolyData* highlight_poly_ = NULL;

    void BuildLargeGlyph();
    void BuildSmallGlyph();
};

#endif