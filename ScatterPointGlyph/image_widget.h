/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef IMAGE_WIDGET_H_
#define IMAGE_WIDGET_H_


#include <vector>
using namespace std;

#include <vtk3DWidget.h>
#include <vtkCommand.h>
#include <vtkObjectFactory.h>
#include <vtkImagePlaneWidget.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkImageActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellPicker.h>
#include <QObject>
#include <QString>
#include <QColor>
#include <QImage>

class ImageWidget : public vtk3DWidget
{
public:
    static ImageWidget* New();

    vtkTypeMacro(ImageWidget, vtk3DWidget);

    void SetData(int w, int h, vector<float>& rgba_values);

    virtual void SetEnabled(int enabled);
	virtual void PlaceWidget(double bounds[6]) {}

    virtual void Modified();

protected:
    ImageWidget();
    virtual ~ImageWidget();

    int image_width_, image_height_;
    vector<float> image_rgba_values_;
    vtkImageActor* image_actor_ = NULL;

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