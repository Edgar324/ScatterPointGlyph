/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef RADAR_BAND_WIDGET_NO_SALIENCY_H_
#define RADAR_BAND_WIDGET_NO_SALIENCY_H_

#include "glyph_widget.h"

class RadarBandWidgetNoSaliency : public GlyphWidget
{
public:
    static RadarBandWidgetNoSaliency* New();

    vtkTypeMacro(RadarBandWidgetNoSaliency, GlyphWidget);

    virtual void SetEnabled(int enabled);

protected:
    RadarBandWidgetNoSaliency();
    ~RadarBandWidgetNoSaliency();

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