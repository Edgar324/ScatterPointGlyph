/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef DENSITY_RENDERING_WIDGET_H_
#define DENSITY_RENDERING_WIDGET_H_

#include <vtk3DWidget.h>
#include <vector>
using namespace std;

class vtkActor;
class vtkPolyDataMapper;
class vtkPolyData;

class DensityRenderingWidget : public vtk3DWidget
{
public:
    static DensityRenderingWidget* New();

    vtkTypeMacro(DensityRenderingWidget, vtk3DWidget);

    void SetData(vector<vector<float>>& pos, vector<int>& cluster_index);

    virtual void SetEnabled(int enabled);
	virtual void PlaceWidget(double bounds[6]) {}

protected:
    DensityRenderingWidget();
    ~DensityRenderingWidget();

private:
    vtkActor* actor_ = NULL;
    vtkPolyDataMapper* data_mapper_ = NULL;
    vtkPolyData* polydata_ = NULL;

    vector<vector<float>> point_pos_;
    vector<int> cluster_index_;
    vector<float> adaptive_rate_;

    void BuildRepresentation();
};

#endif