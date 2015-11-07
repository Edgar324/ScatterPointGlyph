#ifndef SCATTER_POINT_GLYPH_H
#define SCATTER_POINT_GLYPH_H

#include <QtWidgets/QMainWindow>
#include "ui_scatter_point_glyph.h"

#include <string>
#include <vector>

#include <vtkAutoInit.h>
#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

class QMenu;
class QVTKWidget;
class vtkRenderer;
class vtkUnstructuredGrid;
class vtkActor;

struct RenderingLayer {
	std::string name;
	vtkActor* actor;
};

class ScatterPointGlyph : public QMainWindow
{
	Q_OBJECT

public:
	ScatterPointGlyph(QWidget *parent = 0);
	~ScatterPointGlyph();

protected:
	void contextMenuEvent(QContextMenuEvent *event);

private:
	Ui::ScatterPointGlyphClass ui_;
	QVTKWidget* main_view_;
	vtkRenderer* main_renderer_;
	QMenu* context_menu_;

	vtkUnstructuredGrid* scatter_point_data_;

	std::vector< RenderingLayer > rendering_layer_vec_;

	void InitWidget();
	void AddPointData2View();

private slots:
	void OnActionOpenTriggered();
	void OnActionCloseTriggered();
	void OnActionExitTriggered();
	void OnActionHierarchicalClusteringTriggered();
	void OnActionPerceptionDrivenTriggered();
};

#endif // SCATTER_POINT_GLYPH_H
