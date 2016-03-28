#ifndef SCATTER_POINT_GLYPH_H
#define SCATTER_POINT_GLYPH_H

#include <QtWidgets/QMainWindow>
#include "ui_scatter_point_glyph.h"
#include "ui_map_controller.h"
#include "ui_option_dialog.h"

#include <string>
#include <vector>
#include <QtCore/QTimer>

#include <vtkAutoInit.h>
#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

class QMenu;
class QVTKWidget;
class QDockWidget;
class QActionGroup;
class QSlider;
class QLabel;
class QTableView;
class QStandardItemModel;
class QProgressBar;
class QProgressDialog;
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
	vtkRenderer* indicator_renderer_;
	vtkRenderer* other_renderer_;

	Ui::MapController map_control_ui_;
	QWidget* map_control_widget_;

	ParallelCoordinate* parallel_coordinate_;
	ParallelDataset* parallel_dataset_;
	QDockWidget* parallel_coordinate_panel_;

	PathExploreWidget* path_explore_view_;
	PathDataset* path_dataset;
	QDockWidget* path_explore_panel_;

	TreeMapView* tree_map_view_;
	QDockWidget* tree_map_panel_;

    QStandardItemModel* table_model_;
    QTableView* detailed_data_tableview_;
    QDockWidget* data_table_panel_;

	PointRenderingLayer* original_point_rendering_layer_;
	PointRenderingLayer* un_rendering_layer_;

	TransMap* trans_map_;
	TransMapData* transmap_data_;

	QActionGroup* sys_mode_action_group_;
	QActionGroup* main_view_interaction_mode_group_;
	QActionGroup* transmap_tip_mode_group_;
    QActionGroup* color_mapping_group_;
	QActionGroup* example_action_group_;

	QProgressDialog* progress_dialog_;

	// the only dataset served in the memory
	ScatterPointDataset* scatter_point_dataset_;

	// system status
	SystemMode sys_mode_;
	TreeCommon* cluster_tree_;
	int current_view_level_;
	std::vector< int > var_axis_order;

	// thresholds for each cluster method
	float multi_label_threshold_ = 0.2;
	float ncuts_threshold_ = 0.2;
	int hier_number_threshold_ = 1;

	// This factor is used to indicate how many pixels are used to 
	// map the screen distance to the real distance for level cluster retrieval.
	// This should be adjusted according to different application scenarios
	// according to the complexity of the data.
	// For example, if the data is rather simple, this factor should be small, great otherwise.
	// In our experiment, we use 3 for the UCI dataset and 5 for the meteorological data.
	int label_size_factor_ = 3;

	// size of the glyph in pixel radius
	int glyph_pixel_radius_ = { 50 };

	float time_scale_ = 0.0;
	QTimer move_focus_timer_;

	void InitWidget();
	void AddPointData2View();

	float GetMainViewDisPerPixel();
	void GetSceneRange(float& left, float& right, float& bottom, float& top);
	
	void UpdatePointMap();
	void UpdateParallelCoordinate();
	void UpdateTransmap();
	void UpdatePathMap();
	void UpdateTreemap();
    void UpdateTableView();

	void UpdateMenus();

	void ClearAllViews();
	void UpdateAllViews();

	void MoveMainViewToFocus();

	void LoadScData(QString file_path);
	void LoadGscData(QString file_path);

private slots:
	void OnActionOpenVtkFileTriggered();
	void OnActionOpenScatterFileTriggered();
    void OnActionOpenGscFileTriggered();
	void OnActionCloseTriggered();
	void OnActionExitTriggered();
	void OnActionOpenExampleDataTriggered();

	void OnSysmodeChanged();
	void OnSysOptionTriggered();

	void OnBeginClusteringTriggered();
	void OnExecClusteringTriggered();
	void OnClusterFinished();
	void OnMainViewUpdated();
	void OnViewLevelChanged();

	void OnGlyphSelected(int x, int y);
	void OnMainviewLeftButtonUp();
	void OnMainViewRightButtonDown();
	void OnMouseDragmove(int x, int y);;
	void OnSplitClusterTriggered();
	void OnMergeClusterTriggered();
	void OnTreemapNodeSelected(int node_id);
	void OnTransmapHighlightVarChanged(int var_index);
	void OnPcpHighlightVarChanged(int var_index);
	void OnTreemapHighlightVarChagned(int var_index);

	void OnMainViewInteractionModeChanged();
	void OnSavePathSequenceTriggered();
	void OnShowMstTriggered();
	void OnShowVarTrendTriggered();
    void OnMappingVarValueTriggered();

    void OnActionShowScatterPlotTriggered();
	void OnActionShowTransmapTriggered();
    void OnActionShowDensityMapTriggered();
	void OnActionShowTreemapTriggered();
	void OnActionShowTableLensTriggerd();
	void OnActionShowParallelCoordinateTriggered();
    void OnActionShowDataTableTriggered();
    void OnActionShowMapTriggered();

    void OnActionShowResult2DSPTriggered();
    void OnActionShowResult3DSpTriggered();

	void OnGlyphSizeChanged();
	void OnFocusTimerRunOut();
};

#endif // SCATTER_POINT_GLYPH_H
