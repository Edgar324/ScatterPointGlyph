/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef INDICATOR_WIDGET_H_
#define INDICATOR_WIDGET_H_

#include <vtk3DWidget.h>
#include <vtkTextActor3D.h>
#include <vector>
using namespace std;

class vtkActor;
class vtkPolyDataMapper;
class vtkPolyData;

class IndicatorWidget : public vtk3DWidget
{
public:
    static IndicatorWidget* New();

    vtkTypeMacro(IndicatorWidget, vtk3DWidget);

    void SetData(int max_count);
    void SetRenderer(vtkRenderer* renderer);

    virtual void SetEnabled(int enabled);
	virtual void PlaceWidget(double bounds[6]) {}

protected:
    IndicatorWidget();
    ~IndicatorWidget();

private:
    int max_count_ = 0;

    vtkRenderer* renderer_ = nullptr;

    vtkActor* actor_ = NULL;
    vtkPolyDataMapper* data_mapper_ = NULL;
    vtkPolyData* poly_ = NULL;

    vtkTextActor3D* indicator_text_;

    void BuildRepresentation();
};

#endif