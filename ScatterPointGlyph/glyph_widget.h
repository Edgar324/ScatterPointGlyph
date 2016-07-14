/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef GLYPH_WIDGET_H_
#define GLYPH_WIDGET_H_

#include <vector>
using namespace std;

#include <vtk3DWidget.h>
#include <vtkCommand.h>
#include <vtkObjectFactory.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellPicker.h>
#include <QObject>
#include <QString>
#include <QColor>
#include "glyph_dataset.h"
#include "glyph_rendering_widget.h"

class GlyphWidget : public vtk3DWidget
{
public:
    static GlyphWidget* New();

    vtkTypeMacro(GlyphWidget, vtk3DWidget);

    void SetData(GlyphObject* object);
    void SetRenderingWidget(GlyphRenderingWidget* widget);

    virtual void SetEnabled(int enabled);
	virtual void PlaceWidget(double bounds[6]) {}

    virtual void Modified();

protected:
    GlyphWidget();
    virtual ~GlyphWidget();

    enum GlyphState {
        NORMAL = 0x0,
        SELECTION
    };

    const int SEG_PER_PIE = 10;
    const int SEG_PER_CIRCLE = 100;
    const float PIE_VAL = 3.14159265; 
    GlyphObject* glyph_object_ = NULL;
    GlyphRenderingWidget* rendering_widget_ = NULL;
    
    GlyphState glyph_state_ = NORMAL;

    vtkCellPicker* picker_ = NULL;

    virtual void BuildRepresentation();
    virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();

    static void ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata);    
};

#endif