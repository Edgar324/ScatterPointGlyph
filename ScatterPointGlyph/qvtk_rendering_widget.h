#ifndef QVTK_RENDERING_WIDGET_H_
#define QVTK_RENDERING_WIDGET_H_

#include <vtkAutoInit.h>
#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <vector>
using namespace std;

#include <QVTKWidget.h>
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

class QTimer;
class CNode;
class vtkRenderer;

class QVtkRenderingWidget : public QVTKWidget
{
	Q_OBJECT

public:
	QVtkRenderingWidget();
	~QVtkRenderingWidget();

	void ShowTooltip(QString tip);
	void HideTooltip();
    void SetBrushEnabled(bool enabled);

    void GetViewport(float& left, float& right, float& bottom, float& top);
    float GetDistancePerPixel();
    void GetBrushingPath(vector<float>& path) { path = brushing_path; }
    vtkRenderer* renderer() { return main_renderer_; }

signals:
    void ViewUpdated();
    void SelectionChanged();

protected:
	void wheelEvent(QWheelEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
	void mousePressEvent(QMouseEvent* event);
	void mouseReleaseEvent(QMouseEvent* event);

private:
	QTimer* timer_;
	bool is_wheel_updated_;

	QPoint global_mouse_pos_;

    vtkRenderer* main_renderer_;
    //vtkRenderer* indicator_renderer_;

    bool is_brush_enabled_ = false;
    vector<float> brushing_path;

    vtkActor* selection_brush_actor;
	vtkPolyDataMapper* selection_brush_mapper;
	vtkPolyData* selection_brush_poly;

	private slots:
		void OnTimeout();
};

#endif