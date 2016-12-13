/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef POINT_RENDERING_WIDGET_H_
#define POINT_RENDERING_WIDGET_H_

#include <vector>
using namespace std;
#include <vtk3DWidget.h>
#include <QtGui/QColor>

class MultivariateDataset;
class vtkActor;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkScalarBarActor;
class vtkPVScalarBarActor;
class vtkLookupTable;
class vtkRenderer;

class PointRenderingWidget : public vtk3DWidget
{
public:
    static PointRenderingWidget* New();

	vtkTypeMacro(PointRenderingWidget, vtk3DWidget);
	void PrintSelf(ostream& os, vtkIndent indent) {}

	virtual void SetEnabled(int);
	virtual void PlaceWidget(double bounds[6]) {}

	void SetData(MultivariateDataset* data);
    void SetViewpot(float left, float right, float bottom, float top);
    void SetSelectedIds(vector<vector<int>>& ids);
    void SetColorMappingOn(QString name, vector<float>& vals, float min_val, float max_val);
    void SetColorMappingOff();

protected:
    PointRenderingWidget();
	~PointRenderingWidget();

private:
	MultivariateDataset* dataset_ = NULL;
    vector<vector<int>> selected_ids_;

    float view_left_ = 0, view_right_ = 1, view_bottom_ = 0, view_top_ = 1;

	vtkActor* actor_;
	vtkPolyDataMapper* mapper_;
	vtkPolyData* poly_data_;

    vtkPVScalarBarActor* bar_actor_;
    vtkLookupTable* scalar_lookup_table_;

    void BuildPointRepresentation();
    void BuildDensityRepresentation();
};

#endif