#ifndef SCATTER_POINT_GLYPH_H
#define SCATTER_POINT_GLYPH_H

#include <QtWidgets/QMainWindow>
#include "ui_scatter_point_glyph.h"
#include "ui_map_controller.h"
#include "ui_option_dialog.h"
#include "ui_projection_control.h"

#include <string>
#include <vector>
using namespace std;
#include <QtCore/QTimer>

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

class ScatterPointDataset;
class GlyphRenderingWidget;
class GlyphDataset;
class ParallelCoordinate;
class ParallelDataset;

class PointRenderingLayer;

class QVtkRenderingWidget;
class TreeCommon;

class TableLens;

class VolumeRenderer;

class TransMap;
class TransMapData;
class PathExploreWidget; 
class PathDataset;
class TreeMap;
class VariableSelectionWidget;

class ScatterPointGlyph : public QMainWindow
{
	Q_OBJECT

public:
	ScatterPointGlyph(QWidget *parent = 0);
	~ScatterPointGlyph();

	enum ClusteringMode {
		HIER_MODE = 0x0,
		CHAMELEON_MODE,
		NCUTS_MODE,
		MULTI_LABEL_MODE,
        VIEW_DEPENDENT_MODE,
        CLUSTER_PROJECTION_MODE
	};

protected:


private:
	Ui::ScatterPointGlyphClass ui_;
    Ui::ProjectionControlWidget ui_projection_control_;

    // main view for glyph rendering
    GlyphRenderingWidget* glyph_widget_;
    GlyphDataset* glyph_dataset_;

    // Projection control widget
    QWidget* projection_control_widget_;

    // parallel coordinate widget
    ParallelCoordinate* parallel_coordinate_;
	ParallelDataset* parallel_dataset_;
	QDockWidget* parallel_coordinate_panel_;

    // Tree map view
    TreeMap* tree_map_;
	QDockWidget* tree_map_panel_;

    // Table lens view
    TableLens* table_lens_;
    QDockWidget* table_lens_panel_;

    // Data table view
    QStandardItemModel* detailed_data_model_;
    QTableView* detailed_data_tableview_;
    QDockWidget* data_table_panel_;

    // Variable selection widget
    VariableSelectionWidget* variable_selection_widget_;
    QDockWidget* var_selection_panel_;

    // Volume rendering widget
    VolumeRenderer* volume_renderer_;
    vector<float> voxel_values_;
    QDockWidget* volume_render_panel_;

    // the only dataset served in the memory
	ScatterPointDataset* scatter_point_dataset_ = NULL;
    vector<int> selected_var_index_;

    QActionGroup* example_action_group_;
    QActionGroup* color_mapping_group_;
    QActionGroup* clustering_mode_action_group_;
    QActionGroup* glyph_view_interaction_mode_group_;

    ClusteringMode clustering_mode_ = VIEW_DEPENDENT_MODE;
    // thresholds for each cluster method
	float multi_label_threshold_ = 0.1;
	float ncuts_threshold_ = 0.1;
	int hier_number_threshold_ = 1;
    float glyph_size_ = 50;

	TreeCommon* cluster_tree_ = NULL;
    vector<int> selected_cluster_ids_;
    vector<vector<int>> selected_point_ids_;

    // Setup new widgets for the views
	void InitWidget();
    // Initialize views for exploration after a dataset is loaded
    void InitExploration();
    void InitAllViews();
    void UpdateMenus();

    void UpdateAllViews();
    void ClearAllViews();

    void UpdateGlyphWidget();
    void UpdatePcp();
    void UpdateTreemap();
    void UpdateTableLens();
    void UpdateDataTable();
    void UpdateVolumeRender();

	/*vtkRenderer* main_renderer_;
	vtkRenderer* indicator_renderer_;
	vtkRenderer* other_renderer_;

	Ui::MapController map_control_ui_;
	QWidget* map_control_widget_;

	PathExploreWidget* path_explore_view_;
	PathDataset* path_dataset;
	QDockWidget* path_explore_panel_;

	PointRenderingLayer* original_point_rendering_layer_;
	PointRenderingLayer* un_rendering_layer_;

	TransMap* trans_map_;
	TransMapData* transmap_data_;*/

	
	//QActionGroup* transmap_tip_mode_group_;
	

	//QProgressDialog* progress_dialog_;

	//// system status
	//int current_view_level_;
	//std::vector<int> var_axis_order;

	//// This factor is used to indicate how many pixels are used to 
	//// map the screen distance to the real distance for level cluster retrieval.
	//// This should be adjusted according to different application scenarios
	//// according to the complexity of the data.
	//// For example, if the data is rather simple, this factor should be small, great otherwise.
	//// In our experiment, we use 3 for the UCI dataset and 5 for the meteorological data.
	//int label_size_factor_ = 3;

	//// size of the glyph in pixel radius
	//int glyph_pixel_radius_ = { 50 };

	//float time_scale_ = 0.0;
	//QTimer move_focus_timer_;


	/*void AddPointData2View();

	float GetMainViewDisPerPixel();
	
	void UpdatePointMap();
	
	void UpdateTransmap();
	void UpdatePathMap();

	void MoveMainViewToFocus();*/

private slots:
	void OnActionOpenVtkFileTriggered();
	void OnActionOpenScatterFileTriggered();
    void OnActionOpenGscFileTriggered();
	void OnActionCloseTriggered();
	void OnActionExitTriggered();
	void OnActionOpenExampleDataTriggered();

    void OnApplyProjectionTriggered();
    void OnProjectionMethodChanged();

	void OnClusteringModeChanged();
    void OnClusteringOptionTriggered();

    void OnGlyphViewportUpdated();
    void OnGlyphSelectionChanged();
    void OnBrushSelectionChanged();
    void OnVariableSelectionChanged();

    void OnSplitClusterTriggered();
    void OnRecursiveSplittingTriggered();
	void OnMergeClusterTriggered();

    void OnVariableMappingTriggered();
    void OnGlyphViewInteractionModeChanged();

    void OnActionShowScatterPlotTriggered();
	void OnActionShowGlyphTriggered();
    void OnActionShowDensityMapTriggered();
	void OnActionShowTreemapTriggered();
	void OnActionShowTableLensTriggerd();
	void OnActionShowParallelCoordinateTriggered();
    void OnActionShowVolumeRendererTriggered();
    void OnActionShowMapTriggered();
    void OnActionShowDataTableTriggered();

    void OnActionEvaluateQualityTriggered();

    void OnActionShowMeanTriggered();
    void OnActionShowVarianceTriggered();


	/*void OnBeginClusteringTriggered();
	void OnClusterFinished();
	void OnViewLevelChanged();

	void OnGlyphSelected(int x, int y);
	void OnMainviewLeftButtonUp();
	void OnMainViewRightButtonDown();
	void OnMouseDragmove(int x, int y);
	void OnTreemapNodeSelected(int node_id);
	void OnTransmapHighlightVarChanged(int var_index);
	void OnPcpHighlightVarChanged(int var_index);
	void OnTreemapHighlightVarChagned(int var_index);

	void OnSavePathSequenceTriggered();
	void OnShowMstTriggered();
	void OnShowVarTrendTriggered();

    void OnActionShowResult2DSPTriggered();
    void OnActionShowResult3DSpTriggered();

	void OnGlyphSizeChanged();
	void OnFocusTimerRunOut();
    */
};

#endif // SCATTER_POINT_GLYPH_H
