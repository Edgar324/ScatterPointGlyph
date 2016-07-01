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
#include <QObject>
#include <QString>

class GlyphWidget : public QObject, public vtk3DWidget
{
    Q_OBJECT

public:
    static GlyphWidget* New();

    vtkTypeMacro(GlyphWidget, vtk3DWidget);

    void SetData(int cluster_index, vector<QString>& names, 
        vector<float>& means, vector<float>& std_devs, 
        float saliency, int node_count, int max_node_count);

    virtual void SetEnabled(int);
	virtual void PlaceWidget(double bounds[6]) {}
	void PlaceWidget() {
		this->Superclass::PlaceWidget();
	}
	void PlaceWidget(double xmin, double xmax, double ymin, double ymax,
		double zmin, double zmax) {
		this->Superclass::PlaceWidget(xmin, xmax, ymin, ymax, zmin, zmax);
	}

protected:
    GlyphWidget();
    virtual ~GlyphWidget();

    virtual void BuildRepresentation();
    virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();

private:
    int cluster_index_ = -1;
    int node_count_ = 0; 
    int max_node_count_ = 1;
    float saliency_ = 0;
    vector<QString> names_;
    vector<float> means_;
    vector<float> std_devs_;

    static void ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata);    
};

#endif