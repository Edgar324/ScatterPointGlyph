/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef MAP_RENDERING_WIDGET_H_
#define MAP_RENDERING_WIDGET_H_

#include <vtk3DWidget.h>
#include <vector>
using namespace std;

class vtkActor;
class vtkPolyDataMapper;
class vtkPolyData;

class MapRenderingWidget : public vtk3DWidget
{
public:
    static MapRenderingWidget* New();

    vtkTypeMacro(MapRenderingWidget, vtk3DWidget);

    virtual void SetEnabled(int enabled);
	virtual void PlaceWidget(double bounds[6]) {}

protected:
    MapRenderingWidget();
    ~MapRenderingWidget();

private:
    vtkActor* actor_ = NULL;
    vtkPolyDataMapper* data_mapper_ = NULL;
    vtkPolyData* polydata_ = NULL;

    bool is_data_loaded_ = false;

    void BuildRepresentation();
};

#endif