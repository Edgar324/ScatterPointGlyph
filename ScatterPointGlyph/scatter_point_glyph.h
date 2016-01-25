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
class QActionGroup;
class vtkRenderer;
class vtkUnstructuredGrid;
class vtkActor;

class PointRenderingLayer;
class ScatterPointDataset;
class ScatterPointView;
class TreeCommon;
class ParallelCoordinate;
class ParallelDataset;
class TransMap;
class TransMapData;
class PathExploreWidget; 
class PathDataset;
class TreeMapView;

class ScatterPointGlyph : public QMainWindow
{
	Q_OBJECT

public:
	ScatterPointGlyph(QWidget *parent = 0);
	~ScatterPointGlyph();

	enum SystemMode {
		HIER_MODE = 0x0,
		CHAMELEON_MODE,
		NCUTS_MODE,
		MULTI_LABEL_MODE
	};

protected:


private:
	Ui::ScatterPointGlyphClass ui_;
	ScatterPointView* main_view_;
	vtkRenderer* main_renderer_;
	vtkRenderer* dis_matrix_renderer_;
	vtkRenderer* other_renderer_;

	ParallelCoordinate* parallel_coordinate_;
	ParallelDataset* parallel_dataset_;
	QDockWidget* parallel_coordinate_panel_;

	PathExploreWidget* path_explore_view_;
	PathDataset* path_dataset;
	QDockWidget* path_explore_panel_;

	TreeMapView* tree_map_view_;
	QDockWidget* tree_map_panel_;

	PointRenderingLayer* original_point_rendering_layer_;
	PointRenderingLayer* un_rendering_layer_;

	TransMap* trans_map_;
	TransMapData* transmap_data_;

	QActionGroup* sys_mode_action_group_;
	QActionGroup* main_view_interaction_mode_group_;
	QActionGroup* transmap_tip_mode_group_;

	SystemMode sys_mode_;
	TreeCommon* cluster_tree_;

	ScatterPointDataset* scatter_point_dataset_;

	int cluster_num;
	std::vector< int > cluster_index;

	void InitWidget();
	void AddPointData2View();

	float GetMainViewDisPerPixel();
	void GetSceneRange(float& left, float& right, float& bottom, float& top);
	
	void GenerateParallelDataset(ParallelDataset* pdata, std::vector< int >& cluster_ids);
	void UpdateParallelCoordinate();
	void UpdateTransmap();
	void UpdatePathMap();
	void UpdateTreemap();

	void UpdateMenus();

private slots:
	void OnActionOpenVtkFileTriggered();
	void OnActionOpenScatterFileTriggered();
	void OnActionCloseTriggered();
	void OnActionExitTriggered();

	void OnSysmodeChanged();
	void OnExecClusteringTriggered();
	void OnClusterFinished();
	void OnMainViewUpdated();

	void OnGlyphSelected(int x, int y);
	void OnMainviewLeftButtonUp();
	void OnMainViewRightButtonDown();
	void OnMouseDragmove(int x, int y);;
	void OnSplitClusterTriggered();
	void OnMergeClusterTriggered();
	void OnTreemapNodeSelected(int node_id);

	void OnMainViewInteractionModeChanged();
	void OnSavePathSequenceTriggered();
	void OnShowMstTriggered();
	void OnShowVarTrendTriggered();
};

#endif // SCATTER_POINT_GLYPH_H
