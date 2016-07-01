/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef RADAR_BAND_WIDGET_H_
#define RADAR_BAND_WIDGET_H_

#include "glyph_widget.h"

class RadarBandWidget : public GlyphWidget
{
    Q_OBJECT

public:
    static RadarBandWidget* New();

    vtkTypeMacro(RadarBandWidget, GlyphWidget);

protected:
    RadarBandWidget();
    ~RadarBandWidget();

    virtual void BuildRepresentation();
    virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();
};

#endif