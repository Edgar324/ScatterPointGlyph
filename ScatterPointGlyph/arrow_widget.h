/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef ARROW_WIDGET_H_
#define ARROW_WIDGET_H_

#include <vtk3DWidget.h>
#include <vector>
using namespace std;

class vtkActor;
class vtkPolyDataMapper;
class vtkPolyData;

class ArrowWidget : public vtk3DWidget
{
public:
    static ArrowWidget* New();

    vtkTypeMacro(ArrowWidget, vtk3DWidget);

    void SetData(vector<float>& points, float band_size);

    virtual void SetEnabled(int enabled);
	virtual void PlaceWidget(double bounds[6]) {}

protected:
    ArrowWidget();
    ~ArrowWidget();

private:
    vector<float> point_pos_;
    float node_radius_ = 0.1;

    vtkActor* arrow_actor_ = NULL;
    vtkPolyDataMapper* arrow_data_mapper_ = NULL;
    vtkPolyData* arrow_poly_ = NULL;

    void BuildRepresentation();
};

#endif