#include "scatter_point_glyph.h"
#include <fstream>
#include <QVTKWidget.h>
#include <vtkSmartPointer.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkDataObject.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkWriter.h>
#include <vtkCamera.h>
#include <vtkMatrix4x4.h>
#include <vtkInteractorObserver.h>
#include <vtkDelaunay2D.h>

#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QActionGroup>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QSlider>
#include <QtWidgets/QTableView>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QProgressDialog>
#include <QtGui/QStandardItemModel>
#include <QtCore/QTimer>

#include "point_rendering_layer.h"
#include "scatter_point_dataset.h"
#include "scatter_grid_dataset.h"
#include "scatter_point_view.h"
#include "hierarchical_tree.h"
#include "multi_label_tree.h"
#include "ncut_tree.h"
#include "parallel_coordinate.h"
#include "transmap.h"
#include "transmap_data.h"
#include "path_explore_widget.h"
#include "path_dataset.h"
#include "tour_path_generator.h"
#include "tree_map_view.h"
#include "variable_selection_dialog.h"
#include "utility.h"
#include "quality_metric.h"
#include "wrf_data_manager.h"
#include "glyph_design_dialog.h"
#include "var_selection_widget.h"

#define USE_QUALITY_METRIC
//#define SAVE_PROJECTION
//#define USE_SAVED_PROJECTION

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), scatter_point_dataset_(NULL), sys_mode_(MULTI_LABEL_MODE), cluster_tree_(NULL),
	current_view_level_(0) {

	ui_.setupUi(this);

	this->InitWidget();
}

ScatterPointGlyph::~ScatterPointGlyph() {

}

void ScatterPointGlyph::InitWidget() {
	this->setDockOptions(QMainWindow::AllowNestedDocks);

	main_view_ = new ScatterPointView;
	main_view_->setFocusPolicy(Qt::StrongFocus);

	map_control_widget_ = new QWidget;
	map_control_ui_.setupUi(map_control_widget_);
	connect(map_control_ui_.level_slider, SIGNAL(valueChanged(int)), this, SLOT(OnViewLevelChanged()));
	connect(map_control_ui_.glyph_size_spinbox, SIGNAL(valueChanged(int)), this, SLOT(OnGlyphSizeChanged()));

	parallel_coordinate_ = new ParallelCoordinate;
	parallel_coordinate_->setMinimumHeight(200);
	parallel_dataset_ = new ParallelDataset;

	parallel_coordinate_panel_ = new QDockWidget(QString("Parallel Coordinate"), this);
	parallel_coordinate_panel_->setWidget(parallel_coordinate_);
	this->addDockWidget(Qt::BottomDockWidgetArea, parallel_coordinate_panel_);
	parallel_coordinate_panel_->setVisible(false);

    var_selection_widget_ = new VarSelectionWidget;
    var_selection_widget_->setFixedWidth(700);
    var_selection_panel_ = new QDockWidget(QString("Variable Selection"), this);
	var_selection_panel_->setWidget(var_selection_widget_);
	this->addDockWidget(Qt::RightDockWidgetArea, var_selection_panel_);
	var_selection_panel_->setVisible(true);

	tree_map_view_ = new TreeMapView;
	tree_map_view_->setMinimumWidth(700);
	tree_map_panel_ = new QDockWidget(QString("Tree and Table Lens"), this);
	tree_map_panel_->setWidget(tree_map_view_);
	this->tabifyDockWidget(var_selection_panel_, tree_map_panel_);
	tree_map_panel_->setVisible(false);

	path_explore_view_ = new PathExploreWidget;
	path_explore_view_->setMinimumWidth(700);
	path_dataset = new PathDataset;

	path_explore_panel_ = new QDockWidget(QString("Path Records"), this);
	path_explore_panel_->setWidget(path_explore_view_);
	this->tabifyDockWidget(var_selection_panel_, path_explore_panel_);
	path_explore_panel_->setVisible(false);

    detailed_data_tableview_ = new QTableView;
    detailed_data_tableview_->setMinimumWidth(700);
    table_model_ = new QStandardItemModel;
    data_table_panel_ = new QDockWidget(QString("Data"), this);
    data_table_panel_->setWidget(detailed_data_tableview_);
    this->tabifyDockWidget(var_selection_panel_, data_table_panel_);
    data_table_panel_->setVisible(false);

	QVBoxLayout* main_layout = new QVBoxLayout;
	main_layout->addWidget(main_view_);
	main_layout->addWidget(map_control_widget_);
	ui_.centralWidget->setLayout(main_layout);

	main_renderer_ = vtkRenderer::New();
	main_view_->GetRenderWindow()->AddRenderer(main_renderer_);
	main_renderer_->SetViewport(0.0, 0.0, 1.0, 1.0);
	main_renderer_->SetBackground(1.0, 1.0, 1.0);

    indicator_renderer_ = vtkRenderer::New();
    main_view_->GetRenderWindow()->AddRenderer(indicator_renderer_);
	indicator_renderer_->SetViewport(0.8, 0.8, 1.0, 1.0);
	indicator_renderer_->SetBackground(1.0, 1.0, 1.0);

    other_renderer_ = vtkRenderer::New();
    main_view_->GetRenderWindow()->AddRenderer(other_renderer_);
	other_renderer_->SetViewport(1.0, 1.0, 1.0, 1.0);
	other_renderer_->SetBackground(1.0, 1.0, 1.0);

	connect(main_view_, SIGNAL(ViewUpdated()), this, SLOT(OnMainViewUpdated()));
	connect(main_view_, SIGNAL(GlyphSelected(int, int)), this, SLOT(OnGlyphSelected(int, int)));
	connect(main_view_, SIGNAL(LeftButtonUp()), this, SLOT(OnMainviewLeftButtonUp()));
	connect(main_view_, SIGNAL(MouseDrag(int, int)), this, SLOT(OnMouseDragmove(int, int)));
	connect(main_view_, SIGNAL(HighlightVarChanged(int)), this, SLOT(OnTransmapHighlightVarChanged(int)));
	connect(parallel_coordinate_, SIGNAL(HighlightVarChanged(int)), this, SLOT(OnPcpHighlightVarChanged(int)));
	connect(tree_map_view_, SIGNAL(HighlightVarChanged(int)), this, SLOT(OnTreemapHighlightVarChagned(int)));
    connect(tree_map_view_, SIGNAL(NodeSelected(int)), this, SLOT(OnTreemapNodeSelected(int)));
    connect(var_selection_widget_, SIGNAL(SelectionChanged()), this, SLOT(OnVarSelectionChanged()));

	trans_map_ = new TransMap(main_view_);
	trans_map_->SetInteractor(main_view_->GetInteractor());
	trans_map_->SetDefaultRenderer(this->main_renderer_);
    trans_map_->SetIndicatorRenderer(this->indicator_renderer_);
	transmap_data_ = new TransMapData;

	original_point_rendering_layer_ = new PointRenderingLayer;
	original_point_rendering_layer_->SetInteractor(main_view_->GetInteractor());
	original_point_rendering_layer_->SetDefaultRenderer(this->main_renderer_);
    original_point_rendering_layer_->SetColorBarRenderer(this->other_renderer_);
	original_point_rendering_layer_->SetEnabled(true);

	progress_dialog_ = new QProgressDialog;
	progress_dialog_->setWindowModality(Qt::WindowModal);
	progress_dialog_->setLabelText("Processing");
	progress_dialog_->setCancelButton(0);
	progress_dialog_->setMinimumDuration(0);
	progress_dialog_->setRange(0, 0);

	move_focus_timer_.setSingleShot(true);
	connect(&move_focus_timer_, SIGNAL(timeout()), this, SLOT(OnFocusTimerRunOut()));

	sys_mode_action_group_ = new QActionGroup(this);
	sys_mode_action_group_->addAction(ui_.action_hierarchical_clustering);
	sys_mode_action_group_->addAction(ui_.actionChameleon_Clustering);
	sys_mode_action_group_->addAction(ui_.actionNCuts);
	sys_mode_action_group_->addAction(ui_.actionMulti_Label);
	sys_mode_action_group_->setExclusive(true);
	ui_.actionMulti_Label->setChecked(true);
	connect(sys_mode_action_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnSysmodeChanged()));

	main_view_interaction_mode_group_ = new QActionGroup(this);
	main_view_interaction_mode_group_->addAction(ui_.actionSingle_Selection);
	main_view_interaction_mode_group_->addAction(ui_.actionSequence_Selection);
	main_view_interaction_mode_group_->addAction(ui_.actionSelect_Minimum_Path);
	main_view_interaction_mode_group_->addAction(ui_.actionBrush_Path_Sequence);
	main_view_interaction_mode_group_->addAction(ui_.actionBrush_Cluster);
	main_view_interaction_mode_group_->setExclusive(true);
	connect(main_view_interaction_mode_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnMainViewInteractionModeChanged()));

	transmap_tip_mode_group_ = new QActionGroup(this);
	transmap_tip_mode_group_->addAction(ui_.actionOff);
	transmap_tip_mode_group_->setExclusive(true);
	connect(transmap_tip_mode_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnShowVarTrendTriggered()));

    color_mapping_group_ = new QActionGroup(this);
	color_mapping_group_->addAction(ui_.actionColor_Mapping_Off);
	color_mapping_group_->setExclusive(true);
	connect(color_mapping_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnMappingVarValueTriggered()));

	example_action_group_ = new QActionGroup(this);
	example_action_group_->addAction(ui_.actionIris);
	example_action_group_->addAction(ui_.actionWine);
	example_action_group_->addAction(ui_.actionAuto_MPG);
	example_action_group_->addAction(ui_.actionWdbc);
	example_action_group_->addAction(ui_.actionMeteo_Case);
	example_action_group_->setExclusive(true);
	connect(example_action_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnActionOpenExampleDataTriggered()));

	// load data actions
    //connect(ui_.actionOpen_File, SIGNAL(triggered()), this, SLOT(OnActionOpenGridFileTriggered()));
    connect(ui_.actionOpen_File, SIGNAL(triggered()), this, SLOT(OnActionOpenScatterFileTriggered()));
	connect(ui_.action_close, SIGNAL(triggered()), this, SLOT(OnActionCloseTriggered()));
	connect(ui_.actionExit, SIGNAL(triggered()), this, SLOT(OnActionExitTriggered()));
    

	// actions for tips on the cluster transition map
    ui_.mainToolBar->insertAction(ui_.actionShow_Minimum_Spanning_Tree, ui_.menuColor_Mapping->menuAction());
	ui_.mainToolBar->insertAction(ui_.actionShow_Minimum_Spanning_Tree, ui_.menuShow_Sequence->menuAction());
	connect(ui_.actionShow_Minimum_Spanning_Tree, SIGNAL(triggered()), this, SLOT(OnShowMstTriggered()));

	// actions for manipulating cluster nodes
	connect(ui_.actionExec, SIGNAL(triggered()), this, SLOT(OnExecClusteringTriggered()));
	connect(ui_.actionBegin_Clustering, SIGNAL(triggered()), this, SLOT(OnBeginClusteringTriggered()));
	connect(ui_.actionSplit, SIGNAL(triggered()), this, SLOT(OnSplitClusterTriggered()));
	connect(ui_.actionMerge, SIGNAL(triggered()), this, SLOT(OnMergeClusterTriggered()));

	// action for saving exploration path
	connect(ui_.actionAdd_Path_Sequence, SIGNAL(triggered()), this, SLOT(OnSavePathSequenceTriggered()));

	// action for view visibility
    connect(ui_.actionShow_Scatter_Plot, SIGNAL(triggered()), this, SLOT(OnActionShowScatterPlotTriggered()));
	connect(ui_.actionShow_Transmap, SIGNAL(triggered()), this, SLOT(OnActionShowTransmapTriggered()));
	connect(ui_.actionShow_Tree_Map, SIGNAL(triggered()), this, SLOT(OnActionShowTreemapTriggered()));
	connect(ui_.actionShow_Table_Lens, SIGNAL(triggered()), this, SLOT(OnActionShowTableLensTriggerd()));
	connect(ui_.actionShow_PCP, SIGNAL(triggered()), this, SLOT(OnActionShowParallelCoordinateTriggered()));
    connect(ui_.actionShow_Density_Map, SIGNAL(triggered()), this, SLOT(OnActionShowDensityMapTriggered()));
    connect(ui_.actionShow_Data_Table, SIGNAL(triggered()), this, SLOT(OnActionShowDataTableTriggered()));
    connect(ui_.actionShow_Map, SIGNAL(triggered()), this, SLOT(OnActionShowMapTriggered()));

    connect(ui_.actionShow_2D_Result_Scatter_Plot, SIGNAL(triggered()), this, SLOT(OnActionShowResult2DSPTriggered()));
    connect(ui_.actionShow_Result_3D_Scatter_Plot, SIGNAL(triggered()), this, SLOT(OnActionShowResult3DSpTriggered()));

	connect(ui_.actionClusterOptions, SIGNAL(triggered()), this, SLOT(OnSysOptionTriggered()));
}

void ScatterPointGlyph::OnActionOpenVtkFileTriggered() {
	if (scatter_point_dataset_ == NULL) scatter_point_dataset_ = new ScatterPointDataset;
	scatter_point_dataset_->ClearData();

	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName("./TestData/timestep0_0003600.vtk");
	reader->SetReadAllScalars(1);
	reader->SetReadAllVectors(1);
	reader->Update();
	vtkDataObject* data = reader->GetOutput();
	vtkUnstructuredGrid* scatter_point_data = vtkUnstructuredGrid::SafeDownCast(data);

	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer< vtkPolyData >::New();
	polydata->SetPoints(scatter_point_data->GetPoints());

	vtkIntArray* id_array = vtkIntArray::SafeDownCast(scatter_point_data->GetPointData()->GetArray(0));
	vtkFloatArray* speed_array = vtkFloatArray::SafeDownCast(scatter_point_data->GetPointData()->GetArray(1));
	vtkFloatArray* xparray = vtkFloatArray::SafeDownCast(scatter_point_data->GetPointData()->GetArray(2));
	vtkFloatArray* yparray = vtkFloatArray::SafeDownCast(scatter_point_data->GetPointData()->GetArray(3));

	float x_scale1 = 0.35, x_scale2 = 0.48;
	float y_scale1 = 0.65, y_scale2 = 0.77;
	double bounds[6];
	scatter_point_data->GetBounds(bounds);

	for (int i = 0; i < scatter_point_data->GetPoints()->GetNumberOfPoints(); ++i) {
		float x = scatter_point_data->GetPoint(i)[0];
		float y = scatter_point_data->GetPoint(i)[1];
		if (x > bounds[0] + (bounds[1] - bounds[0]) * x_scale1 && x < bounds[0] + (bounds[1] - bounds[0]) * x_scale2
			&& y > bounds[2] + (bounds[3] - bounds[2]) * y_scale1 && y < bounds[2] + (bounds[3] - bounds[2]) * y_scale2) {
			std::vector<float> pos;
			std::vector<float> value;
			pos.push_back(x);
			pos.push_back(y);
			value.push_back(speed_array->GetValue(i));
			value.push_back(xparray->GetValue(i));
			value.push_back(yparray->GetValue(i));

			scatter_point_dataset_->original_point_pos.push_back(pos);
			scatter_point_dataset_->original_point_values.push_back(value);
		}
	}

	scatter_point_dataset_->DirectConstruct();
	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionOpenScatterFileTriggered() {
	QString file_path = QFileDialog::getOpenFileName(this, tr("Open Scatter Point Data"), ".", "*.sc");
	if (file_path.length() == 0) return;

	LoadScData(file_path);

	this->UpdateMenus();

	this->AddPointData2View();

    var_selection_widget_->SetData(scatter_point_dataset_->var_names, scatter_point_dataset_->var_colors);
    selected_var_index_.resize(scatter_point_dataset_->var_num);
    for (int i = 0; i < scatter_point_dataset_->var_num; ++i) selected_var_index_[i] = i;

	this->label_size_factor_ = 3.0;
}

void ScatterPointGlyph::OnActionOpenGscFileTriggered() {
	QString file_path = QFileDialog::getOpenFileName(this, tr("Open Scatter Point Data"), ".", "*.gsc");
	if (file_path.length() == 0) return;

	LoadGscData(file_path);

	this->UpdateMenus();

	this->AddPointData2View();

    var_selection_widget_->SetData(scatter_point_dataset_->var_names, scatter_point_dataset_->var_colors);
    selected_var_index_.resize(scatter_point_dataset_->var_num);
    for (int i = 0; i < scatter_point_dataset_->var_num; ++i) selected_var_index_[i] = i;

	this->label_size_factor_ = 5.0;
}

void ScatterPointGlyph::OnActionOpenExampleDataTriggered() {
	this->OnActionCloseTriggered();

	QString file_path = "./TestData/";
	if (ui_.actionIris->isChecked())
		file_path += "iris.sc";
	else if (ui_.actionAuto_MPG->isChecked())
		file_path += "auto-mpg.sc";
	else if (ui_.actionWine->isChecked())
		file_path += "wine.sc";
	else if (ui_.actionWdbc->isChecked())
		file_path += "wdbc.sc";
	else if (ui_.actionMeteo_Case->isChecked())
        //file_path += "plot.gsc";
		file_path += "agent_step100_status.gsc";

	if (!ui_.actionMeteo_Case->isChecked()) {
		LoadScData(file_path);

		QString mds_path = file_path.left(file_path.length() - 2);
		mds_path = mds_path + QString("mds");
		std::ifstream input(mds_path.toLocal8Bit().data());
		if (input.good()) {
			for (int i = 0; i < scatter_point_dataset_->original_point_pos.size(); ++i)
				for (int j = 0; j < scatter_point_dataset_->original_point_pos[i].size(); ++j)
					input >> scatter_point_dataset_->original_point_pos[i][j];
		}
		input.close();
		
		scatter_point_dataset_->DirectConstruct();

		this->UpdateMenus();

		this->AddPointData2View();

		this->label_size_factor_ = 3.0;
	} else {
		LoadGscData(file_path);

		this->UpdateMenus();

		this->AddPointData2View();

		this->label_size_factor_ = 5.0;
	}

    var_selection_widget_->SetData(scatter_point_dataset_->var_names, scatter_point_dataset_->var_colors);
    selected_var_index_.resize(scatter_point_dataset_->var_num);
    for (int i = 0; i < scatter_point_dataset_->var_num; ++i) selected_var_index_[i] = i;
}

void ScatterPointGlyph::LoadScData(QString file_path) {
	if (scatter_point_dataset_ == NULL) scatter_point_dataset_ = new ScatterPointDataset;
	scatter_point_dataset_->ClearData();

	//std::ifstream input_file(file_path.toLocal8Bit());
	//std::ifstream input_file("./TestData/auto-mpg.sc");
	//std::ifstream input_file("./TestData/wine.sc");
	//std::ifstream input_file("./TestData/iris.sc");
	//std::ifstream input_file("./TestData/wdbc.sc");
	//std::ifstream input_file("./TestData/yeast.sc");
	std::ifstream input_file(file_path.toLocal8Bit().data());
	char char_str[1000];
	input_file.getline(char_str, 1000);
	QString value_str = QString::fromLocal8Bit(char_str);
	QStringList value_list = value_str.split(' ');

	int record_num, var_num;
	record_num = value_list.at(0).toInt();
	var_num = value_list.at(1).toInt();

	input_file.getline(char_str, 1000);
	value_str = QString::fromLocal8Bit(char_str);
	value_list = value_str.split(' ');
	for (int i = 0; i < value_list.size(); ++i) scatter_point_dataset_->var_names.push_back(value_list.at(i));

	VariableSelectionDialog var_dialog;
	var_dialog.SetDatasetInfo(record_num, var_num, scatter_point_dataset_->var_names);
	if (var_dialog.exec() != QDialog::Accepted) {
		input_file.close();
		return;
	}

	scatter_point_dataset_->original_point_pos.resize(record_num);
	scatter_point_dataset_->original_point_values.resize(record_num);
	for (int i = 0; i < record_num; ++i) {
		scatter_point_dataset_->original_point_values[i].resize(var_num);
		scatter_point_dataset_->original_point_pos[i].resize(2);

		input_file.getline(char_str, 1000);
		value_str = QString::fromLocal8Bit(char_str);
		value_list = value_str.split(',');

		for (int j = 0; j < var_num; ++j)
			scatter_point_dataset_->original_point_values[i][j] = value_list.at(j).toFloat();
	}
	input_file.close();

	if (var_dialog.IsAutomaticDimReduction()) {
		int dim_num = var_dialog.GetDimNumber();
		if (dim_num <= 0 && dim_num > var_num) return;
		scatter_point_dataset_->AutoDimReduction(dim_num);
	}
	else {
		std::vector<float> dim_weights;
		std::vector<bool> is_dim_selected;
		var_dialog.GetSelectionResult(dim_weights);
		scatter_point_dataset_->var_weights.clear();
		for (int i = 0; i < dim_weights.size(); ++i)
			if (dim_weights[i] >= 0) {
				is_dim_selected.push_back(true);
			}
			else
				is_dim_selected.push_back(false);
		scatter_point_dataset_->var_weights = dim_weights;
		scatter_point_dataset_->ManualSelectDim(is_dim_selected);
	}

	scatter_point_dataset_->ExecMds();

#ifdef SAVE_PROJECTION
	std::ofstream output("./TestData/temp.mds");
	if (output.good()) {
		for (int i = 0; i < scatter_point_dataset_->original_point_pos.size(); ++i)
			for (int j = 0; j < scatter_point_dataset_->original_point_pos[i].size(); ++j)
				output << scatter_point_dataset_->original_point_pos[i][j] << " ";
	}
	output.close();
#endif

#ifdef USE_SAVED_PROJECTION
	QString mds_path = file_path.left(file_path.length() - 2);
	mds_path = mds_path + QString("mds");
	std::ifstream input(mds_path.toLocal8Bit().data());
	if (input.good()) {
		for (int i = 0; i < scatter_point_dataset_->original_point_pos.size(); ++i)
			for (int j = 0; j < scatter_point_dataset_->original_point_pos[i].size(); ++j)
				input >> scatter_point_dataset_->original_point_pos[i][j];
	}
	input.close();
#endif

	scatter_point_dataset_->DirectConstruct();
}

void ScatterPointGlyph::LoadGscData(QString file_path) {
	if (scatter_point_dataset_ == NULL) scatter_point_dataset_ = new ScatterPointDataset;
	scatter_point_dataset_->ClearData();

	//std::ifstream input_file("./TestData/plot.gsc");
	//std::ifstream input_file("./TestData/veri.gsc");
	std::ifstream input_file(file_path.toLocal8Bit().data());
	char char_str[1000];
	input_file.getline(char_str, 1000);
	QString value_str = QString::fromLocal8Bit(char_str);
	QStringList value_list = value_str.split(' ');

	int record_num, var_num;
	record_num = value_list.at(0).toInt();
	var_num = value_list.at(1).toInt();

	input_file.getline(char_str, 1000);
	value_str = QString::fromLocal8Bit(char_str);
	value_list = value_str.split(' ');
	for (int i = 0; i < value_list.size(); ++i) scatter_point_dataset_->var_names.push_back(value_list.at(i));

	VariableSelectionDialog var_dialog;
	var_dialog.SetDatasetInfo(record_num, var_num, scatter_point_dataset_->var_names);
	if (var_dialog.exec() != QDialog::Accepted) {
		input_file.close();
		return;
	}

	scatter_point_dataset_->original_point_pos.resize(record_num);
	scatter_point_dataset_->original_point_values.resize(record_num);
	for (int i = 0; i < record_num; ++i) {
		scatter_point_dataset_->original_point_values[i].resize(var_num);
		scatter_point_dataset_->original_point_pos[i].resize(2);

		input_file.getline(char_str, 1000);
		value_str = QString::fromLocal8Bit(char_str);
		value_list = value_str.split(',');

		for (int j = 0; j < var_num; ++j)
			scatter_point_dataset_->original_point_values[i][j] = value_list.at(2 + j).toFloat();

		scatter_point_dataset_->original_point_pos[i][0] = value_list.at(0).toFloat();
		scatter_point_dataset_->original_point_pos[i][1] = value_list.at(1).toFloat();
	}
	input_file.close();

	if (var_dialog.IsAutomaticDimReduction()) {
		int dim_num = var_dialog.GetDimNumber();
		if (dim_num <= 0 && dim_num > var_num) return;
		scatter_point_dataset_->AutoDimReduction(dim_num);
	}
	else {
		std::vector<float> dim_weights;
		std::vector<bool> is_dim_selected;
		var_dialog.GetSelectionResult(dim_weights);
		scatter_point_dataset_->var_weights.clear();
		for (int i = 0; i < dim_weights.size(); ++i)
			if (dim_weights[i] >= 0) {
				scatter_point_dataset_->var_weights.push_back(dim_weights[i]);
				is_dim_selected.push_back(true);
			}
			else
				is_dim_selected.push_back(false);
		scatter_point_dataset_->ManualSelectDim(is_dim_selected);
	}

	scatter_point_dataset_->DirectConstruct();
}

void ScatterPointGlyph::OnActionCloseTriggered() {
	this->ClearAllViews();

	this->original_point_rendering_layer_->ClearView();

	if (scatter_point_dataset_ != NULL) delete scatter_point_dataset_;
	scatter_point_dataset_ = NULL;

	if (cluster_tree_ != NULL ) delete cluster_tree_;
	cluster_tree_ = NULL;

	this->main_view_->update();
}

void ScatterPointGlyph::OnActionExitTriggered() {
	// TODO: clear all the memory

	exit(0);
}

void ScatterPointGlyph::OnSysmodeChanged() {
	SystemMode current_mode = sys_mode_;

	if (cluster_tree_ != NULL) {
		QMessageBox msgbox;
		msgbox.setText("The tree structure will be destroyed if you change the clustering mode!");
		msgbox.setInformativeText("Do you want to continue?");
		msgbox.setStandardButtons(QMessageBox::Yes | QMessageBox::Cancel);
		int ret = msgbox.exec();
		switch (ret)
		{
		case QMessageBox::Yes:
			if (ui_.action_hierarchical_clustering->isChecked())
				current_mode = HIER_MODE;
			else if (ui_.actionChameleon_Clustering->isChecked())
				current_mode = CHAMELEON_MODE;
			else if (ui_.actionNCuts->isChecked())
				current_mode = NCUTS_MODE;
			else if (ui_.actionMulti_Label)
				current_mode = MULTI_LABEL_MODE;
			break;
		default:
			switch (current_mode)
			{
			case ScatterPointGlyph::HIER_MODE:
				ui_.action_hierarchical_clustering->setChecked(true);
				break;
			case ScatterPointGlyph::CHAMELEON_MODE:
				ui_.action_hierarchical_clustering->setChecked(true);
				break;
			case ScatterPointGlyph::NCUTS_MODE:
				ui_.actionNCuts->setChecked(true);
				break;
			case ScatterPointGlyph::MULTI_LABEL_MODE:
				ui_.actionMulti_Label->setChecked(true);
				break;
			default:
				break;
			}
			return;			
		}
	} else {
		if (ui_.action_hierarchical_clustering->isChecked())
			current_mode = HIER_MODE;
		else if (ui_.actionChameleon_Clustering->isChecked())
			current_mode = CHAMELEON_MODE;
		else if (ui_.actionNCuts->isChecked())
			current_mode = NCUTS_MODE;
		else if (ui_.actionMulti_Label)
			current_mode = MULTI_LABEL_MODE;
	}

	if (current_mode != sys_mode_) {
		if (cluster_tree_ != NULL) {
			this->ClearAllViews();
			delete cluster_tree_;
			cluster_tree_ = NULL;
		}
		sys_mode_ = current_mode;
	}
}

void ScatterPointGlyph::OnSysOptionTriggered() {
	QDialog option_dialog;
	Ui::OptionDialog option_ui;
	option_ui.setupUi(&option_dialog);
	option_ui.ml_threshold_spinbox->setValue(multi_label_threshold_);
	option_ui.ncuts_threshold_spinbox->setValue(ncuts_threshold_);
	option_ui.hier_number_spinbox->setValue(hier_number_threshold_);
	option_ui.factor_spinbox->setValue(label_size_factor_);

	if (option_dialog.exec() == QDialog::Accepted) {
		multi_label_threshold_ = option_ui.ml_threshold_spinbox->value();
		ncuts_threshold_ = option_ui.ncuts_threshold_spinbox->value();
		hier_number_threshold_ = option_ui.hier_number_spinbox->value();
		label_size_factor_ = option_ui.factor_spinbox->value();
	}
}

void ScatterPointGlyph::OnExecClusteringTriggered() {
	if (scatter_point_dataset_ == NULL) {
		QMessageBox::information(this, tr("Warning"), tr("Please load data first."));
		return;
	}

	ui_.statusBar->showMessage("Execute clustering, progressing ...");
	progress_dialog_->show();

	float dis_per_pixel = this->GetMainViewDisPerPixel();
	
	switch (sys_mode_)
	{
	case ScatterPointGlyph::HIER_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new HierarchicalTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
			connect(cluster_tree_, SIGNAL(finished()), progress_dialog_, SLOT(cancel()));
		}

		HierarchicalTree* hier_tree = dynamic_cast<HierarchicalTree*>(cluster_tree_);
        hier_tree->SetTreeMode(TreeCommon::VIEWING_MODE);
		if (hier_tree != NULL) {
			//hier_tree->SetExpectedClusterNum(1);
			hier_tree->start();
		}
	}
		break;
	case ScatterPointGlyph::CHAMELEON_MODE:
	{
		
	}
		break;
	case ScatterPointGlyph::NCUTS_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new NCutTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
			connect(cluster_tree_, SIGNAL(finished()), progress_dialog_, SLOT(cancel()));
		}

		NCutTree* ncut_tree = dynamic_cast<NCutTree*>(cluster_tree_);
		ncut_tree->SetTreeMode(TreeCommon::VIEWING_MODE);
		if (ncut_tree != NULL) {
			ncut_tree->start();
		}
	}
		break;
	case ScatterPointGlyph::MULTI_LABEL_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new MultiLabelTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
			connect(cluster_tree_, SIGNAL(finished()), progress_dialog_, SLOT(cancel()));
		}

		MultiLabelTree* un_tree = dynamic_cast< MultiLabelTree* >(cluster_tree_);
		un_tree->SetTreeMode(TreeCommon::VIEWING_MODE);
		if (un_tree != NULL) {
			float dis_per_pixel = this->GetMainViewDisPerPixel();
			//un_tree->SetRadiusThreshold(200.0 * dis_per_pixel / (scatter_point_dataset_->original_pos_ranges[0][1] - scatter_point_dataset_->original_pos_ranges[0][0]));
			un_tree->start();
		}
	}
		break;
	default:
		break;
	}
}

void ScatterPointGlyph::OnBeginClusteringTriggered() {
	if (scatter_point_dataset_ == NULL) {
		QMessageBox::information(this, tr("Warning"), tr("Please load data first."));
		return;
	}

	ui_.statusBar->showMessage("Begin clustering, progressing ...");
	progress_dialog_->show();

	float dis_per_pixel = this->GetMainViewDisPerPixel();
	
	switch (sys_mode_)
	{
	case ScatterPointGlyph::NCUTS_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new NCutTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
			connect(cluster_tree_, SIGNAL(finished()), progress_dialog_, SLOT(cancel()));
		}

		NCutTree* ncut_tree = dynamic_cast<NCutTree*>(cluster_tree_);
		ncut_tree->SetTreeMode(TreeCommon::EXPLORATION_MODE);
		if (ncut_tree != NULL) {
			ncut_tree->start();
		}
	}
		break;
	case ScatterPointGlyph::MULTI_LABEL_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new MultiLabelTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
			connect(cluster_tree_, SIGNAL(finished()), progress_dialog_, SLOT(cancel()));
		}

		MultiLabelTree* un_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);
		un_tree->SetTreeMode(TreeCommon::EXPLORATION_MODE);
		if (un_tree != NULL) un_tree->start();
	}
		break;
	default:
		break;
	}
}

void ScatterPointGlyph::AddPointData2View() {
	original_point_rendering_layer_->SetData(scatter_point_dataset_);
	original_point_rendering_layer_->SetEnabled(true);

	main_renderer_->ResetCamera();

	this->UpdateParallelCoordinate();

	main_view_->update();
}

void ScatterPointGlyph::OnClusterFinished() {
	ui_.statusBar->showMessage("Cluster finished.");

	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* multi_label_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);
		if (multi_label_tree == NULL) return;

		float dis_per_pixel = this->GetMainViewDisPerPixel();
		current_view_level_ = multi_label_tree->GetRadiusLevel(glyph_pixel_radius_ * label_size_factor_ * dis_per_pixel / scatter_point_dataset_->max_pos_range);
		map_control_ui_.level_slider->setValue(current_view_level_);
		map_control_ui_.level_index_label->setText(QString("%0").arg(current_view_level_));
	} else {
		current_view_level_ = 1;
	}

	int max_level = cluster_tree_->GetMaxLevel();
    if (current_view_level_ >= max_level) current_view_level_ = max_level - 1;
	map_control_ui_.level_slider->setRange(0, max_level);
	map_control_ui_.level_slider->setValue(current_view_level_);
	map_control_ui_.level_index_label->setText(QString("%0").arg(current_view_level_));
	map_control_ui_.level_name_label->setText(QString("Level(0~%0): ").arg(max_level));

#ifdef USE_QUALITY_METRIC
	QualityMetric* metric = new QualityMetric;
	metric->GenerateQualityMeasures(cluster_tree_);
	metric->SaveMeasures("quality.txt");
	delete metric;
#endif

    current_view_level_ = 0;
	this->UpdateAllViews();
}

void ScatterPointGlyph::OnMainViewUpdated() {
    if (scatter_point_dataset_ == NULL || cluster_tree_ == NULL) return;

	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* multi_label_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);
		if (multi_label_tree == NULL) return;

		float dis_per_pixel = this->GetMainViewDisPerPixel();
		current_view_level_ = multi_label_tree->GetRadiusLevel(glyph_pixel_radius_ * label_size_factor_ * dis_per_pixel / scatter_point_dataset_->max_pos_range);
        int max_level = multi_label_tree->GetMaxLevel();
        if (current_view_level_ >= max_level) current_view_level_ = max_level - 1;
		map_control_ui_.level_slider->setValue(current_view_level_);
		map_control_ui_.level_index_label->setText(QString("%0").arg(current_view_level_));
	}

	this->UpdateAllViews();
}

void ScatterPointGlyph::OnViewLevelChanged() {
	current_view_level_ = map_control_ui_.level_slider->value();
	map_control_ui_.level_index_label->setText(QString("%0").arg(map_control_ui_.level_slider->value()));

	this->UpdateAllViews();
}

void ScatterPointGlyph::UpdateAllViews() {
	if (scatter_point_dataset_ == NULL || cluster_tree_ == NULL) return;
	this->UpdateTransmap();
	this->UpdateParallelCoordinate();
	this->UpdateTreemap();
	this->UpdatePointMap();
    this->UpdateTableView();
}

void ScatterPointGlyph::ClearAllViews() {
	trans_map_->ClearView();
	transmap_data_->ClearData();

	parallel_dataset_->ClearData();
	this->parallel_coordinate_->SetDataset(parallel_dataset_);

	this->tree_map_view_->ClearView();
}

void ScatterPointGlyph::UpdateTableView() {
    std::vector<int> selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);

    std::vector<int> cluster_index;
	cluster_index.resize(scatter_point_dataset_->point_num, -1);
    for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i) {
		std::vector<int> temp_vec;
		cluster_tree_->Traverse(transmap_data_->cluster_nodes[i], temp_vec);
		for (int j = 0; j < temp_vec.size(); ++j) cluster_index[temp_vec[j]] = i;
	}

    table_model_->clear();
    QStringList headers;
    for (int i = 0; i < scatter_point_dataset_->var_num; ++i)
        headers << scatter_point_dataset_->var_names[i];
    table_model_->setHorizontalHeaderLabels(headers);

    for (int i = 0; i < selection_index.size(); ++i) {
		for (int j = 0; j < scatter_point_dataset_->point_num; ++j)
			if (cluster_index[j] == selection_index[i]) {
                int row_count = table_model_->rowCount();
                table_model_->insertRow(row_count);
                for (int k = 0; k < scatter_point_dataset_->var_num; ++k) {
                    table_model_->setData(table_model_->index(row_count, k), scatter_point_dataset_->original_point_values[j][k], Qt::DisplayRole);
                }
			}
	}
    detailed_data_tableview_->setModel(table_model_);
    detailed_data_tableview_->setSortingEnabled(true);
    detailed_data_tableview_->update();
}

void ScatterPointGlyph::UpdateParallelCoordinate() {
	if (transmap_data_ == NULL || cluster_tree_ == NULL || scatter_point_dataset_ == NULL) return;

	if (transmap_data_->cluster_nodes.size() == 0) {
		parallel_dataset_->ClearData();
		parallel_dataset_->subset_names.push_back(QString("MDS Data"));
		parallel_dataset_->subset_records.resize(1);
        std::vector<std::vector<float>> value_ranges;
        value_ranges.resize(selected_var_index_.size());
        for (int i = 0; i < selected_var_index_.size(); ++i) {
            value_ranges[i].resize(2);
            value_ranges[i][0] = 1e20;
            value_ranges[i][1] = -1e20;
        }

		for (int i = 0; i < scatter_point_dataset_->original_point_values.size(); ++i) {
			ParallelRecord* record = new ParallelRecord;
            for (int j = 0; j < selected_var_index_.size(); ++j)
                record->values.push_back(scatter_point_dataset_->original_point_values[i][selected_var_index_[j]]);
            for (int j = 0; j < selected_var_index_.size(); ++j) {
                if (value_ranges[j][0] > record->values[j]) value_ranges[j][0] = record->values[j];
                if (value_ranges[j][1] < record->values[j]) value_ranges[j][1] = record->values[j];
            }
			parallel_dataset_->subset_records[0].push_back(record);
		}

        for (int i = 0; i < selected_var_index_.size(); ++i)
            if (value_ranges[i][1] - value_ranges[i][0] == 0) {
                value_ranges[i][0] -= 0.5;
                value_ranges[i][1] += 0.5;
            }
        for (int i = 0; i < parallel_dataset_->subset_records[0].size(); ++i) {
            for (int j = 0; j < selected_var_index_.size(); ++j) {
                parallel_dataset_->subset_records[0][i]->values[j] = (parallel_dataset_->subset_records[0][i]->values[j] - value_ranges[j][0]) / (value_ranges[j][1] - value_ranges[j][0]);
            }
        }

        parallel_dataset_->axis_anchors.resize(scatter_point_dataset_->var_num);
		for (int i = 0; i < selected_var_index_.size(); ++i) {
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(value_ranges[i][0]));
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(value_ranges[i][1]));
			parallel_dataset_->axis_names.push_back(scatter_point_dataset_->var_names[selected_var_index_[i]]);
		}

		parallel_dataset_->CompleteInput();
		parallel_dataset_->UpdateGaussian();
		parallel_coordinate_->SetDataset(parallel_dataset_);
		Utility::GenerateAxisOrder(parallel_dataset_, var_axis_order);
		parallel_coordinate_->SetAxisOrder(var_axis_order);

		QString str = QString("%0 points.").arg(scatter_point_dataset_->original_point_pos.size());
		ui_.statusBar->showMessage(str);
	} else {
		std::vector<int> selected_cluster_index;
		trans_map_->GetSelectedClusterIndex(selected_cluster_index);
		/*if (selection_index.size() == 0) {
			for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i) selection_index.push_back(i);
		}*/

		// fix the dataset in case of rendering during the data update
		parallel_dataset_->is_updating = true;
		parallel_dataset_->ClearData();

		std::vector<int> cluster_color = original_point_rendering_layer_->GetClusterColor();
		parallel_dataset_->subset_names.resize(selected_cluster_index.size());
		parallel_dataset_->subset_colors.resize(selected_cluster_index.size());
		parallel_dataset_->subset_records.resize(selected_cluster_index.size());
		parallel_dataset_->axis_names.resize(selected_var_index_.size());
		parallel_dataset_->axis_anchors.resize(selected_var_index_.size());

		int cluster_num;
		std::vector<int> cluster_index;
		cluster_index.resize(scatter_point_dataset_->point_num, -1);
		for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i) {
			std::vector<int> temp_vec;
			cluster_tree_->Traverse(transmap_data_->cluster_nodes[i], temp_vec);
			for (int j = 0; j < temp_vec.size(); ++j) cluster_index[temp_vec[j]] = i;
		}
		//cluster_tree_->GetClusterResult(current_view_level_, cluster_num, cluster_index);

        std::vector<std::vector<float>> value_ranges;
        value_ranges.resize(selected_var_index_.size());
        for (int i = 0; i < selected_var_index_.size(); ++i) {
            value_ranges[i].resize(2);
            value_ranges[i][0] = 1e20;
            value_ranges[i][1] = -1e20;
        }

		for (int i = 0; i < selected_cluster_index.size(); ++i) {
			parallel_dataset_->subset_names[i] = QString("Cluster %0").arg(i);
			parallel_dataset_->subset_colors[i] = QColor(cluster_color[3 * selected_cluster_index[i]], cluster_color[3 * selected_cluster_index[i] + 1], cluster_color[3 * selected_cluster_index[i] + 2]);
			for (int j = 0; j < scatter_point_dataset_->point_num; ++j)
				if (cluster_index[j] == selected_cluster_index[i]) {
					ParallelRecord* record = new ParallelRecord;
                    for (int k = 0; k < selected_var_index_.size(); ++k)
                        record->values.push_back(scatter_point_dataset_->original_point_values[j][selected_var_index_[k]]);

                    for (int k = 0; k < record->values.size(); ++k) {
                        if (value_ranges[k][0] > record->values[k]) value_ranges[k][0] = record->values[k];
                        if (value_ranges[k][1] < record->values[k]) value_ranges[k][1] = record->values[k];
                    }

					parallel_dataset_->subset_records[i].push_back(record);
				}
		}

        for (int i = 0; i < selected_var_index_.size(); ++i)
            if (value_ranges[i][1] - value_ranges[i][0] == 0) {
                value_ranges[i][0] -= 0.5;
                value_ranges[i][1] += 0.5;
            }
        for (int i = 0; i < parallel_dataset_->subset_records.size(); ++i) {
            for (int j = 0; j < parallel_dataset_->subset_records[i].size(); ++j) {
                for (int k = 0; k < selected_var_index_.size(); ++k) {
                    parallel_dataset_->subset_records[i][j]->values[k] = (parallel_dataset_->subset_records[i][j]->values[k] - value_ranges[k][0]) / (value_ranges[k][1] - value_ranges[k][0]);
                }
            }
        }

        if (selected_cluster_index.size() == 0) {
            for (int i = 0; i < selected_var_index_.size(); ++i) {
                value_ranges[i].resize(2);
                value_ranges[i][0] = scatter_point_dataset_->original_value_ranges[selected_var_index_[i]][0];
                value_ranges[i][1] = scatter_point_dataset_->original_value_ranges[selected_var_index_[i]][1];
            }
        }
		for (int i = 0; i < selected_var_index_.size(); ++i) {
			parallel_dataset_->axis_names[i] = scatter_point_dataset_->var_names[selected_var_index_[i]];
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(value_ranges[i][0]));
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(value_ranges[i][1]));
		}
		parallel_dataset_->CompleteInput();
		parallel_dataset_->UpdateGaussian();

		parallel_dataset_->is_updating = false;

		parallel_coordinate_->SetDataset(parallel_dataset_);
		Utility::GenerateAxisOrder(parallel_dataset_, var_axis_order);

        /*std::vector<int> focus_index;
        if (trans_map_ != NULL) trans_map_->GetFocusVarIndex(focus_index);
        for (int i = 0; i < focus_index.size(); ++i) {
            int temp_pos = -1;
            for (int j = 0; j < var_axis_order.size(); ++j)
                if (var_axis_order[j] == focus_index[i]) {
                    temp_pos = j;
                    break;
                }
            if (temp_pos != -1) {
                for (int j = temp_pos; j > i; --j) {
                    var_axis_order[j] = var_axis_order[j - 1];
                }
                var_axis_order[i] = focus_index[i];
            }
        }*/

		parallel_coordinate_->SetAxisOrder(var_axis_order);
		parallel_coordinate_->update();

		int selected_count = 0;
		for (int i = 0; i < parallel_dataset_->subset_records.size(); ++i)
			selected_count += parallel_dataset_->subset_records[i].size();
		QString str = QString("%0 out of %1 (%2%) points are selected.").arg(selected_count).arg(scatter_point_dataset_->original_point_pos.size()).arg((int)selected_count * 100 / scatter_point_dataset_->original_point_pos.size());
		ui_.statusBar->showMessage(str);
	}
}

void ScatterPointGlyph::UpdateTransmap() {
	if (transmap_data_ == NULL || cluster_tree_ == NULL || scatter_point_dataset_ == NULL) return;

	for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
		transmap_data_->cluster_nodes[i]->is_highlighted = false;
	transmap_data_->ClearData();

	int cluster_num;
	std::vector<int> cluster_index;
	float dis_per_pixel = this->GetMainViewDisPerPixel();
	// dis_per_pixel * 100.0
	cluster_tree_->GetClusterResult(current_view_level_, transmap_data_->cluster_nodes);
	cluster_tree_->GetClusterResult(current_view_level_, cluster_num, cluster_index);
	
	std::vector<QColor > colors;
	for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
		colors.push_back(transmap_data_->cluster_nodes[i]->color);

	original_point_rendering_layer_->SetClusterIndex(cluster_num, cluster_index, colors);
	original_point_rendering_layer_->SetHighlightCluster(-1);

	transmap_data_->cluster_num = cluster_num;
	transmap_data_->var_num = scatter_point_dataset_->var_weights.size();
	transmap_data_->dataset = scatter_point_dataset_;
	transmap_data_->ProcessData();

	trans_map_->SetNodeRadius(dis_per_pixel * glyph_pixel_radius_);
	trans_map_->SetData(scatter_point_dataset_, transmap_data_, selected_var_index_);
	trans_map_->SetAxisOrder(var_axis_order);
    std::vector<QColor > point_colors;
    for (int i = 0; i < cluster_index.size(); ++i) {
        point_colors.push_back(colors[cluster_index[i]]);
    }
    trans_map_->UpdateDensityActor(point_colors);

	if (ui_.actionShow_Transmap->isChecked()) {
		trans_map_->SetEnabled(true);
	}

	main_view_->update();
}


void ScatterPointGlyph::UpdatePointMap() {
	if (trans_map_ == NULL) return;

	std::vector<int> selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}
}

void ScatterPointGlyph::UpdatePathMap() {

}

void ScatterPointGlyph::UpdateTreemap() {
	if (trans_map_ == NULL || cluster_tree_ == NULL) return;

	std::vector<int> selected_ids;
	trans_map_->GetSelectedClusterIds(selected_ids);

	//UncertaintyTree* un_tree = dynamic_cast<UncertaintyTree*>(cluster_tree_);
	cluster_tree_->SortTree(selected_ids);

	std::vector<int> selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	std::vector<bool> is_selected;
	is_selected.resize(transmap_data_->cluster_nodes.size(), false);

	std::vector<CNode*> selected_nodes;
	if (selected_ids.size() != 0) {
        std::vector<QString > names;
        std::vector<QColor > colors;
        for (int i = 0; i < selected_var_index_.size(); ++i) {
            names.push_back(scatter_point_dataset_->var_names[selected_var_index_[i]]);
            colors.push_back(scatter_point_dataset_->var_colors[selected_var_index_[i]]);
        }

		for (int i = 0; i < selection_index.size(); ++i) {
			selected_nodes.push_back(transmap_data_->cluster_nodes[selection_index[i]]);
			is_selected[selection_index[i]] = true;
		}
		for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
			if (!is_selected[i]) selected_nodes.push_back(transmap_data_->cluster_nodes[i]);
		tree_map_view_->SetData(cluster_tree_, selected_var_index_, selected_nodes, selected_ids.size(), var_axis_order, names, colors);
	}
	else {
        std::vector<QString > names;
        std::vector<QColor > colors;
        for (int i = 0; i < selected_var_index_.size(); ++i) {
            names.push_back(scatter_point_dataset_->var_names[selected_var_index_[i]]);
            colors.push_back(scatter_point_dataset_->var_colors[selected_var_index_[i]]);
        }
		selected_nodes = transmap_data_->cluster_nodes;
		tree_map_view_->SetData(cluster_tree_, selected_var_index_, selected_nodes, selected_nodes.size(), var_axis_order, names, colors);
	}

	tree_map_view_->scene()->update();
}

float ScatterPointGlyph::GetMainViewDisPerPixel() {
	double point_one[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 1, 0, 0, point_one);
	double point_two[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 2, 0, 0, point_two);

	return (1.0 / abs(point_one[0] - point_two[0]));
}

void ScatterPointGlyph::GetSceneRange(float& left, float& right, float& bottom, float& top) {
	double viewport[4];
	this->main_renderer_->GetViewport(viewport);
	double point_one[4];
	vtkInteractorObserver::ComputeDisplayToWorld(this->main_renderer_, viewport[0] * this->main_view_->width(), viewport[1] * this->main_view_->height(), 0, point_one);
	double point_two[4];
	vtkInteractorObserver::ComputeDisplayToWorld(this->main_renderer_, viewport[2] * this->main_view_->width(), viewport[3] * this->main_view_->height(), 0, point_two);
	left = point_one[0];
	right = point_two[0];
	bottom = point_one[1];
	top = point_two[1];
}

void ScatterPointGlyph::OnGlyphSelected(int x, int y) {
	if (scatter_point_dataset_ == NULL || cluster_tree_ == NULL) return;

	this->UpdateParallelCoordinate();
    this->UpdateTableView();
	trans_map_->SetAxisOrder(var_axis_order);
	this->UpdateTreemap();
	this->UpdatePointMap();
	this->main_view_->update();
}

void ScatterPointGlyph::OnMainviewLeftButtonUp() {
	if (trans_map_ != NULL) {
		trans_map_->OnMouseReleased();
		if (ui_.actionBrush_Cluster->isChecked() || ui_.actionBrush_Path_Sequence->isChecked()) {
			this->UpdateParallelCoordinate();
			trans_map_->SetAxisOrder(var_axis_order);
			this->UpdateTreemap();
			this->UpdatePointMap();
			this->main_view_->update();
		}
	}
}

void ScatterPointGlyph::OnMainViewRightButtonDown() {

}

void ScatterPointGlyph::OnMouseDragmove(int x, int y) {
	if (trans_map_ != NULL) trans_map_->OnMouseMove(x, y);
}

void ScatterPointGlyph::OnSplitClusterTriggered() {
	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* multi_label_tree = dynamic_cast< MultiLabelTree* >(cluster_tree_);

		int temp_cluster_index = trans_map_->GetSelectedClusterIndex();
		if (temp_cluster_index != -1) {
			progress_dialog_->show();

			multi_label_tree->SplitCluster(transmap_data_->cluster_nodes[temp_cluster_index]->id());

			progress_dialog_->cancel();

			float dis_per_pixel = this->GetMainViewDisPerPixel();
			current_view_level_ = multi_label_tree->GetRadiusLevel(glyph_pixel_radius_ * label_size_factor_ * dis_per_pixel / scatter_point_dataset_->max_pos_range);
			int max_level = multi_label_tree->GetMaxLevel();
			if (current_view_level_ >= max_level) current_view_level_ = max_level;
			map_control_ui_.level_slider->setValue(current_view_level_);
			map_control_ui_.level_index_label->setText(QString("%0").arg(current_view_level_));

			UpdateTransmap();
			UpdateTreemap();
			UpdateParallelCoordinate();
		}
	}

    if (sys_mode_ == NCUTS_MODE && cluster_tree_ != NULL) {
		NCutTree* ncut_tree = dynamic_cast< NCutTree* >(cluster_tree_);

		int temp_cluster_index = trans_map_->GetSelectedClusterIndex();
		if (temp_cluster_index != -1) {
			ncut_tree->SplitCluster(transmap_data_->cluster_nodes[temp_cluster_index]->id());
			UpdateTransmap();
			UpdateTreemap();
			UpdateParallelCoordinate();
		}
	}

    // update the level slider
	if (cluster_tree_ != NULL) {
		int max_level = cluster_tree_->GetMaxLevel();
		map_control_ui_.level_slider->setRange(0, max_level);
		map_control_ui_.level_slider->setValue(current_view_level_);
		map_control_ui_.level_index_label->setText(QString("%0").arg(current_view_level_));
		map_control_ui_.level_name_label->setText(QString("Level(0~%0): ").arg(max_level));
	}
}

void ScatterPointGlyph::OnMergeClusterTriggered() {
	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* un_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);

		std::vector<int> merged_clusters;
		trans_map_->GetSelectedClusterIndex(merged_clusters);

		std::vector<int> cluster_ids;
		for (int i = 0; i < merged_clusters.size(); ++i)
			cluster_ids.push_back(transmap_data_->cluster_nodes[merged_clusters[i]]->id());
		un_tree->MergeClusters(cluster_ids);

		UpdateTransmap();
		UpdateTreemap();
		UpdateParallelCoordinate();
	}
}

void ScatterPointGlyph::OnTreemapNodeSelected(int node_id) {
	trans_map_->OnNodeSelected(node_id);

	std::vector<int> selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}

	this->UpdateParallelCoordinate();

	this->MoveMainViewToFocus();

	this->main_view_->update();
}

void ScatterPointGlyph::MoveMainViewToFocus() {
	time_scale_ = 0.0;

	this->trans_map_->ForceFocusCenter();
	move_focus_timer_.start(100);
}

void ScatterPointGlyph::OnFocusTimerRunOut() {
	if (time_scale_ < 1.0) {
		time_scale_ += 0.2;

		this->trans_map_->MoveViewToFocus(time_scale_);
		move_focus_timer_.start(100);
	}
}

void ScatterPointGlyph::OnMainViewInteractionModeChanged() {
	if (ui_.actionSingle_Selection->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_SINGLE_CLUSTER);
	} else if (ui_.actionSequence_Selection->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_MULTI_CLUSTERS);
	} else if (ui_.actionSelect_Minimum_Path->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_MINIMUM_PATH);
	} else if (ui_.actionBrush_Path_Sequence->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_BRUSHED_PATH_SEQUENCE);
	} else if (ui_.actionBrush_Cluster->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_BRUSHED_CLUSTERS);
	} else {
		trans_map_->SetInteractionState(TransMap::NORMAL);
	}
}

void ScatterPointGlyph::OnSavePathSequenceTriggered() {
	std::list<CNode*> node_seq = trans_map_->GetNodeSequence();

	PathRecord* record = new PathRecord;
	record->item_values.resize(node_seq.size());
	record->item_color.resize(node_seq.size());
	std::list<CNode*>::iterator node_iter = node_seq.begin();

	int count = 0;
	while (node_iter != node_seq.end()) {
		record->item_values[count] = (*node_iter)->average_values;
		record->item_color[count] = QColor(128, 128, 128);
		count++;
		node_iter++;
	}

	record->change_values.resize(node_seq.size() - 1);
	for (int i = 0; i < node_seq.size() - 1; ++i) {
		record->change_values[i].resize(record->item_values[0].size());
		for (int j = 0; j < record->item_values[0].size(); ++j)
			record->change_values[i][j] = record->item_values[i + 1][j] - record->item_values[i][j];
	}
	for (int k = 0; k < record->item_values[0].size(); ++k)
		record->var_names.push_back("Test");

	path_dataset->path_records.push_back(record);
	path_explore_view_->SetData(path_dataset);
}

void ScatterPointGlyph::OnShowMstTriggered() {
	if (ui_.actionShow_Minimum_Spanning_Tree->isChecked()) {
		//ui_.actionShow_Sequence->setChecked(false);
		trans_map_->ShowMinimumSpanningTree(true);
	} else {
		trans_map_->ShowMinimumSpanningTree(false);
	}

	main_view_->update();
}

void ScatterPointGlyph::OnShowVarTrendTriggered() {
	if (!ui_.actionOff->isChecked()) {
		int var_index = -1;
		QList< QAction* > actions = transmap_tip_mode_group_->actions();
		for (int i = 0; i < actions.size(); ++i)
			if (actions.at(i)->isChecked()) {
				var_index = i;
				break;
			}

		trans_map_->ShowVarTrend(var_index - 1);
	} else {
		trans_map_->ShowVarTrend(-1);
	}

	main_view_->update();

	this->UpdateTreemap();
	this->UpdateParallelCoordinate();
}

void ScatterPointGlyph::OnMappingVarValueTriggered() {
    if (!ui_.actionColor_Mapping_Off->isChecked()) {
		int var_index = -1;
		QList< QAction* > actions = color_mapping_group_->actions();
		for (int i = 0; i < actions.size(); ++i)
			if (actions.at(i)->isChecked()) {
				var_index = i;
				break;
			}

		std::vector<float> values;
        for (int i = 0; i < scatter_point_dataset_->original_point_values.size(); ++i) 
            values.push_back(scatter_point_dataset_->original_point_values[i][var_index - 1]);
        original_point_rendering_layer_->SetPointValue(values);
        main_view_->update();
	} else {
		trans_map_->ShowVarTrend(-1);
        std::vector<float> values;
        original_point_rendering_layer_->SetPointValue(values);
        main_view_->update();
	}
}

void ScatterPointGlyph::UpdateMenus() {
	QList< QAction* > actions = transmap_tip_mode_group_->actions(); 
	for (int i = 1; i < actions.size(); ++i) {
		transmap_tip_mode_group_->removeAction(actions.at(i));
		ui_.menuShow_Sequence->removeAction(actions.at(i));
		delete actions.at(i);
	}

	ui_.actionOff->setChecked(true);
	for (int i = 0; i < scatter_point_dataset_->var_names.size(); ++i) {
		QAction* action = ui_.menuShow_Sequence->addAction(scatter_point_dataset_->var_names[i]);
		action->setCheckable(true);
		action->setChecked(false);
		transmap_tip_mode_group_->addAction(action);
	}
	//transmap_tip_mode_group_->setExclusive(true);

    QList< QAction* > actions_mapping = color_mapping_group_->actions(); 
	for (int i = 1; i < actions_mapping.size(); ++i) {
		color_mapping_group_->removeAction(actions_mapping.at(i));
		ui_.menuColor_Mapping->removeAction(actions_mapping.at(i));
		delete actions_mapping.at(i);
	}

	ui_.actionColor_Mapping_Off->setChecked(true);
	for (int i = 0; i < scatter_point_dataset_->var_names.size(); ++i) {
		QAction* action = ui_.menuColor_Mapping->addAction(scatter_point_dataset_->var_names[i]);
		action->setCheckable(true);
		action->setChecked(false);
		color_mapping_group_->addAction(action);
	}
}

void ScatterPointGlyph::OnActionShowScatterPlotTriggered() {
    original_point_rendering_layer_->SetScatterPointEnabled(ui_.actionShow_Scatter_Plot->isChecked());
    this->main_view_->update();
}

void ScatterPointGlyph::OnActionShowTransmapTriggered()
{
	trans_map_->SetEnabled(ui_.actionShow_Transmap->isChecked());
}

void ScatterPointGlyph::OnActionShowTreemapTriggered()
{
	tree_map_view_->SetTreeMapVisible(ui_.actionShow_Tree_Map->isChecked());
	tree_map_panel_->setVisible(tree_map_view_->IsVisible());
}

void ScatterPointGlyph::OnActionShowTableLensTriggerd()
{
	tree_map_view_->SetTableLensVisible(ui_.actionShow_Table_Lens->isChecked());
	tree_map_panel_->setVisible(tree_map_view_->IsVisible());
}

void ScatterPointGlyph::OnActionShowParallelCoordinateTriggered()
{
	parallel_coordinate_panel_->setVisible(ui_.actionShow_PCP->isChecked());
}

void ScatterPointGlyph::OnActionShowDensityMapTriggered() {
    if (trans_map_ != NULL) {
        trans_map_->SetDensityMapVisibility(ui_.actionShow_Density_Map->isChecked());
    }
    if (tree_map_view_ != NULL) {
        tree_map_view_->SetTreeMapUsingColor(ui_.actionShow_Density_Map->isChecked());
    }
}

void ScatterPointGlyph::OnActionShowDataTableTriggered() {
    data_table_panel_->setVisible(ui_.actionShow_Data_Table->isChecked());
}

void ScatterPointGlyph::OnActionShowMapTriggered() {
    if (original_point_rendering_layer_ != NULL) {
        original_point_rendering_layer_->SetMapEnabled(ui_.actionShow_Map->isChecked());
        this->main_view_->update();
    }
}

void ScatterPointGlyph::OnTransmapHighlightVarChanged(int var_index)
{
	parallel_coordinate_->SetHighlightAxis(var_index);
	tree_map_view_->SetHighlightVarIndex(var_index);
}

void ScatterPointGlyph::OnPcpHighlightVarChanged(int var_index)
{
	trans_map_->HighlightVar(var_index);
	tree_map_view_->SetHighlightVarIndex(var_index);
}

void ScatterPointGlyph::OnTreemapHighlightVarChagned(int var_index)
{
	trans_map_->HighlightVar(var_index);
	parallel_coordinate_->SetHighlightAxis(var_index);
}

void ScatterPointGlyph::OnGlyphSizeChanged()
{
	glyph_pixel_radius_ = map_control_ui_.glyph_size_spinbox->value();

	this->UpdateTransmap();
}

void ScatterPointGlyph::OnActionShowResult2DSPTriggered()
{

}

void ScatterPointGlyph::OnActionShowResult3DSpTriggered()
{

}

void ScatterPointGlyph::OnVarSelectionChanged()
{
    vector<bool> is_selected;
    var_selection_widget_->GetSelection(is_selected);
    selected_var_index_.clear();
    for (int i = 0; i < is_selected.size(); ++i)
        if (is_selected[i]) selected_var_index_.push_back(i);
    this->UpdateAllViews();
}
