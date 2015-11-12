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
class QDockWidget;
class vtkRenderer;
class vtkUnstructuredGrid;
class vtkActor;
class LayerControlWidget;
class RenderingLayerModel;
class RenderingLayer;
class HierParaWidget;
class HierSolver;
class ClusterGlyphLayer;
class PointRenderingLayer;
class GestaltProcessor;

class ScatterPointGlyph : public QMainWindow
{
	Q_OBJECT

public:
	ScatterPointGlyph(QWidget *parent = 0);
	~ScatterPointGlyph();

	enum SystemMode {
		HIER_MODE = 0x0,
		PERCEPTION_MODE,
	};

protected:

private:
	Ui::ScatterPointGlyphClass ui_;
	QVTKWidget* main_view_;
	vtkRenderer* main_renderer_;

	QDockWidget* layer_control_panel_;
	LayerControlWidget* layer_control_widget_;
	RenderingLayerModel* rendering_layer_model_;

	QDockWidget* hier_para_panel_;
	HierParaWidget* hier_para_widget_;

	SystemMode sys_mode_;

	HierSolver* hier_solver_;
	int expected_cluster_num_;
	ClusterGlyphLayer* cluster_glyph_layer_;
	PointRenderingLayer* point_rendering_layer_;
	PointRenderingLayer* original_rendering_layer_;

	GestaltProcessor* gestalt_processor_;

	vtkUnstructuredGrid* scatter_point_data_;
	std::vector< std::vector< float > > point_pos_;
	std::vector< std::vector< float > > point_values_;
	std::vector< float > variable_weight_;


	void InitWidget();
	void AddPointData2View();

	void PreProcess();
	void HierarchicalPreProcess();
	void PerceptionPreProcess();

	void ExecPerceptionClustering();

private slots:
	void OnActionOpenTriggered();
	void OnActionCloseTriggered();
	void OnActionExitTriggered();
	void OnActionHierarchicalClusteringTriggered();
	void OnActionPerceptionDrivenTriggered();
	
	void OnHierClusterNumberChanged(int);
	void OnCombinedClusterUpdated(int, int);
	void OnOneStepHierFinished();
	void OnGestaltUpdated();
	void OnKmeansClusterFinished();

	void NormalizePosition(std::vector< std::vector< float > >& vec);
	void NormalizeVector(std::vector< std::vector< float > >& vec);
};

#endif // SCATTER_POINT_GLYPH_H
