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
#include <QtGui/QStandardItemModel>

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

//#define USE_QUALITY_METRIC

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

	tree_map_view_ = new TreeMapView;
	tree_map_view_->setMinimumWidth(700);
	tree_map_panel_ = new QDockWidget(QString("Tree Map and Table Lens"), this);
	tree_map_panel_->setWidget(tree_map_view_);
	this->addDockWidget(Qt::RightDockWidgetArea, tree_map_panel_);
	tree_map_panel_->setVisible(false);

	path_explore_view_ = new PathExploreWidget;
	path_explore_view_->setMinimumWidth(700);
	path_dataset = new PathDataset;

	path_explore_panel_ = new QDockWidget(QString("Path Records"), this);
	path_explore_panel_->setWidget(path_explore_view_);
	this->tabifyDockWidget(tree_map_panel_, path_explore_panel_);
	path_explore_panel_->setVisible(false);

    detailed_data_tableview_ = new QTableView;
    detailed_data_tableview_->setMinimumWidth(700);
    table_model_ = new QStandardItemModel;
    data_table_panel_ = new QDockWidget(QString("Data"), this);
    data_table_panel_->setWidget(detailed_data_tableview_);
    this->tabifyDockWidget(tree_map_panel_, data_table_panel_);
    data_table_panel_->setVisible(false);

	QVBoxLayout* main_layout = new QVBoxLayout;
	main_layout->addWidget(main_view_);
	main_layout->addWidget(map_control_widget_);
	ui_.centralWidget->setLayout(main_layout);

	main_renderer_ = vtkRenderer::New();
	main_view_->GetRenderWindow()->AddRenderer(main_renderer_);
	main_renderer_->SetViewport(0.0, 0.0, 0.85, 1.0);
	main_renderer_->SetBackground(1.0, 1.0, 1.0);

    indicator_renderer_ = vtkRenderer::New();
    main_view_->GetRenderWindow()->AddRenderer(indicator_renderer_);
	indicator_renderer_->SetViewport(0.85, 0.8, 1.0, 1.0);
	indicator_renderer_->SetBackground(1.0, 1.0, 1.0);

    other_renderer_ = vtkRenderer::New();
    main_view_->GetRenderWindow()->AddRenderer(other_renderer_);
	other_renderer_->SetViewport(0.85, 0.0, 1.0, 0.8);
	other_renderer_->SetBackground(1.0, 1.0, 1.0);

	connect(main_view_, SIGNAL(ViewUpdated()), this, SLOT(OnMainViewUpdated()));
	connect(main_view_, SIGNAL(GlyphSelected(int, int)), this, SLOT(OnGlyphSelected(int, int)));
	connect(main_view_, SIGNAL(LeftButtonUp()), this, SLOT(OnMainviewLeftButtonUp()));
	connect(main_view_, SIGNAL(MouseDrag(int, int)), this, SLOT(OnMouseDragmove(int, int)));
	connect(main_view_, SIGNAL(HighlightVarChanged(int)), this, SLOT(OnTransmapHighlightVarChanged(int)));
	connect(parallel_coordinate_, SIGNAL(HighlightVarChanged(int)), this, SLOT(OnPcpHighlightVarChanged(int)));
	connect(tree_map_view_, SIGNAL(HighlightVarChanged(int)), this, SLOT(OnTreemapHighlightVarChagned(int)));
    connect(tree_map_view_, SIGNAL(NodeSelected(int)), this, SLOT(OnTreemapNodeSelected(int)));

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

	// load data actions
    connect(ui_.actionOpen_File, SIGNAL(triggered()), this, SLOT(OnActionOpenGridFileTriggered()));
    //connect(ui_.actionOpen_File, SIGNAL(triggered()), this, SLOT(OnActionOpenScatterFileTriggered()));
    

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
			std::vector< float > pos;
			std::vector< float > value;
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
	//QString file_path = QFileDialog::getOpenFileName(this, tr("Open Scatter Point Data"), ".", "*.sc");
	//if (file_path.length() == 0) return;

	if (scatter_point_dataset_ == NULL) scatter_point_dataset_ = new ScatterPointDataset;
	scatter_point_dataset_->ClearData();

	//std::ifstream input_file(file_path.toLocal8Bit());
	//std::ifstream input_file("./TestData/auto-mpg.sc");
	std::ifstream input_file("./TestData/wine.sc");
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
	} else {
		std::vector< float > dim_weights;
		std::vector< bool > is_dim_selected;
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

	scatter_point_dataset_->ExecMds();
	scatter_point_dataset_->DirectConstruct();

	this->UpdateMenus();

	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionOpenGridFileTriggered() {
    WrfDataManager* manager = WrfDataManager::GetInstance();
    manager->LoadEnsembleData(WRF_ACCUMULATED_PRECIPITATION, std::string("G:/Data/ens_l/apcp_sfc_latlon_all_20150401_20150430_liaoLeeVC4.nc"));
    manager->LoadEnsembleData(WRF_PRECIPITABLE_WATER, std::string("G:/Data/ens_l/pwat_eatm_latlon_all_20150401_20150430_liao5WubTq.nc"));
    manager->LoadEnsembleData(WRF_T2M, std::string("G:/Data/ens_l/tmp_2m_latlon_all_20150401_20150430_liaoSwbPtU.nc"));
    manager->LoadEnsembleData(WRF_MSLP, std::string("G:/Data/ens_l/pres_msl_latlon_all_20150401_20150430_liaoawNRCR.nc"));

    QDateTime temp_datetime = QDateTime(QDate(2015, 4, 5), QTime(0, 0));

    std::vector< WrfGridValueMap* > apcp_maps, prec_maps, t2m_maps, mslp_maps;
    manager->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_ACCUMULATED_PRECIPITATION, 24, apcp_maps);
    manager->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_PRECIPITABLE_WATER, 24, prec_maps);
    manager->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_T2M, 24, t2m_maps);
    manager->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_MSLP, 24, mslp_maps);

    ScatterGridDataset* grid_data = new ScatterGridDataset;
    scatter_point_dataset_ = grid_data;
    MapRange range = apcp_maps[0]->map_range;
    grid_data->w = apcp_maps[0]->map_range.x_grid_number;
    grid_data->h = apcp_maps[0]->map_range.y_grid_number;
    grid_data->point_num = grid_data->w * grid_data->h;
    grid_data->var_num = 4;
    grid_data->var_names.push_back(QString("APCP"));
    grid_data->var_names.push_back(QString("PREW"));
    grid_data->var_names.push_back(QString("T2M"));
    grid_data->var_names.push_back(QString("MSLP"));
    grid_data->original_point_pos.resize(grid_data->point_num);
    grid_data->original_point_values.resize(grid_data->point_num);
    grid_data->var_weights.resize(4, 0.25);

    int accu_index = 0;
    for (int i = 0; i < grid_data->h; ++i)
        for (int j = 0; j < grid_data->w; ++j) {
            grid_data->original_point_pos[accu_index].resize(2);
            grid_data->original_point_pos[accu_index][0] = range.start_x + j * range.x_grid_space;
            grid_data->original_point_pos[accu_index][1] = range.start_y + i * range.y_grid_space;

            grid_data->original_point_values[accu_index].resize(4);
            grid_data->original_point_values[accu_index][0] = apcp_maps[0]->values[accu_index];
            grid_data->original_point_values[accu_index][1] = prec_maps[0]->values[accu_index];
            grid_data->original_point_values[accu_index][2] = t2m_maps[0]->values[accu_index];
            grid_data->original_point_values[accu_index][3] = mslp_maps[0]->values[accu_index];

            accu_index++;
        }

    scatter_point_dataset_->DirectConstruct();

	this->UpdateMenus();

	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionCloseTriggered() {

}

void ScatterPointGlyph::OnActionExitTriggered() {
	
}

void ScatterPointGlyph::OnSysmodeChanged() {
	SystemMode current_mode;
	if (ui_.action_hierarchical_clustering->isChecked())
		current_mode = HIER_MODE;
	else if (ui_.actionChameleon_Clustering->isChecked())
		current_mode = CHAMELEON_MODE;
	else if (ui_.actionNCuts->isChecked())
		current_mode = NCUTS_MODE;
	else if (ui_.actionMulti_Label)
		current_mode = MULTI_LABEL_MODE;

	if (current_mode != sys_mode_) {
		if (cluster_tree_ != NULL) {
			delete cluster_tree_;
			cluster_tree_ = NULL;
		}

		sys_mode_ = current_mode;
	}
}

void ScatterPointGlyph::OnExecClusteringTriggered() {
	if (scatter_point_dataset_ == NULL) {
		QMessageBox::information(this, tr("Warning"), tr("Please load data first."));
		return;
	}

	float dis_per_pixel = this->GetMainViewDisPerPixel();
	
	switch (sys_mode_)
	{
	case ScatterPointGlyph::HIER_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new HierarchicalTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}

		HierarchicalTree* hier_tree = dynamic_cast<HierarchicalTree*>(cluster_tree_);
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

	float dis_per_pixel = this->GetMainViewDisPerPixel();

	switch (sys_mode_)
	{
	case ScatterPointGlyph::NCUTS_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new NCutTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
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
		}

		MultiLabelTree* un_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);
		un_tree->SetTreeMode(TreeCommon::EXPLORATION_MODE);
		if (un_tree != NULL) {
			//float dis_per_pixel = this->GetMainViewDisPerPixel();
			//un_tree->SetRadiusThreshold(200.0 * dis_per_pixel / (scatter_point_dataset_->original_pos_ranges[0][1] - scatter_point_dataset_->original_pos_ranges[0][0]));
			un_tree->start();
		}
	}
	break;
	default:
		break;
	}
}

void ScatterPointGlyph::AddPointData2View() {
	original_point_rendering_layer_->SetData(scatter_point_dataset_);

	main_renderer_->ResetCamera();

	this->UpdateParallelCoordinate();

	main_view_->update();
}

void ScatterPointGlyph::OnClusterFinished() {
	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* multi_label_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);
		if (multi_label_tree == NULL) return;

		float dis_per_pixel = this->GetMainViewDisPerPixel();
		current_view_level_ = multi_label_tree->GetRadiusLevel(label_pixel_radius_ * 3 * dis_per_pixel / scatter_point_dataset_->max_pos_range);
		map_control_ui_.level_slider->setValue(current_view_level_);
		map_control_ui_.level_index_label->setText(QString("%0").arg(current_view_level_));
	} else {
		current_view_level_ = 1;
	}

	int max_level = cluster_tree_->GetMaxLevel();
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

	this->UpdateAllViews();
}

void ScatterPointGlyph::OnMainViewUpdated() {
    if (scatter_point_dataset_ == NULL || cluster_tree_ == NULL) return;

	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* multi_label_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);
		if (multi_label_tree == NULL) return;

		float dis_per_pixel = this->GetMainViewDisPerPixel();
		current_view_level_ = multi_label_tree->GetRadiusLevel(label_pixel_radius_ * 3 * dis_per_pixel / scatter_point_dataset_->max_pos_range);
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
	this->UpdateTransmap();
	this->UpdateParallelCoordinate();
	this->UpdateTreemap();
	this->UpdatePointMap();
    this->UpdateTableView();
}

void ScatterPointGlyph::UpdateTableView() {
    std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);

    std::vector< int > cluster_index;
	cluster_index.resize(scatter_point_dataset_->point_num, -1);
    for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i) {
		std::vector< int > temp_vec;
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
	if (transmap_data_->cluster_nodes.size() == 0) {
		parallel_dataset_->ClearData();
		parallel_dataset_->subset_names.push_back(QString("MDS Data"));
		parallel_dataset_->subset_records.resize(1);
        std::vector< std::vector< float > > value_ranges;
        value_ranges.resize(scatter_point_dataset_->var_num);
        for (int i = 0; i < scatter_point_dataset_->var_num; ++i) {
            value_ranges[i].resize(2);
            value_ranges[i][0] = 1e20;
            value_ranges[i][1] = -1e20;
        }

		for (int i = 0; i < scatter_point_dataset_->original_point_values.size(); ++i) {
			ParallelRecord* record = new ParallelRecord;
			record->values = scatter_point_dataset_->original_point_values[i];
            for (int j = 0; j < record->values.size(); ++j) {
                if (value_ranges[j][0] > record->values[j]) value_ranges[j][0] = record->values[j];
                if (value_ranges[j][1] < record->values[j]) value_ranges[j][1] = record->values[j];
            }
			parallel_dataset_->subset_records[0].push_back(record);
		}

        for (int i = 0; i < scatter_point_dataset_->var_num; ++i)
            if (value_ranges[i][1] - value_ranges[i][0] == 0) {
                value_ranges[i][0] -= 0.5;
                value_ranges[i][1] += 0.5;
            }
        for (int i = 0; i < parallel_dataset_->subset_records[0].size(); ++i) {
            for (int j = 0; j < scatter_point_dataset_->var_num; ++j) {
                parallel_dataset_->subset_records[0][i]->values[j] = (parallel_dataset_->subset_records[0][i]->values[j] - value_ranges[j][0]) / (value_ranges[j][1] - value_ranges[j][0]);
            }
        }

        parallel_dataset_->axis_anchors.resize(scatter_point_dataset_->var_num);
		for (int i = 0; i < scatter_point_dataset_->var_num; ++i) {
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(value_ranges[i][0]));
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(value_ranges[i][1]));
			parallel_dataset_->axis_names.push_back(scatter_point_dataset_->var_names[i]);
		}

		parallel_dataset_->CompleteInput();
		parallel_dataset_->UpdateGaussian();
		parallel_coordinate_->SetDataset(parallel_dataset_);
		Utility::GenerateAxisOrder(parallel_dataset_, var_axis_order);
		parallel_coordinate_->SetAxisOrder(var_axis_order);

		QString str = QString("%0 points.").arg(scatter_point_dataset_->original_point_pos.size());
		ui_.statusBar->showMessage(str);
	} else {
		std::vector< int > selection_index;
		trans_map_->GetSelectedClusterIndex(selection_index);
		/*if (selection_index.size() == 0) {
			for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i) selection_index.push_back(i);
		}*/

		// fix the dataset in case of rendering during the data update
		parallel_dataset_->is_updating = true;
		parallel_dataset_->ClearData();

		std::vector< int > cluster_color = original_point_rendering_layer_->GetClusterColor();
		parallel_dataset_->subset_names.resize(selection_index.size());
		parallel_dataset_->subset_colors.resize(selection_index.size());
		parallel_dataset_->subset_records.resize(selection_index.size());
		parallel_dataset_->axis_names.resize(scatter_point_dataset_->original_value_ranges.size());
		parallel_dataset_->axis_anchors.resize(scatter_point_dataset_->original_value_ranges.size());

		int cluster_num;
		std::vector< int > cluster_index;
		cluster_index.resize(scatter_point_dataset_->point_num, -1);
		for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i) {
			std::vector< int > temp_vec;
			cluster_tree_->Traverse(transmap_data_->cluster_nodes[i], temp_vec);
			for (int j = 0; j < temp_vec.size(); ++j) cluster_index[temp_vec[j]] = i;
		}
		//cluster_tree_->GetClusterResult(current_view_level_, cluster_num, cluster_index);

        std::vector< std::vector< float > > value_ranges;
        value_ranges.resize(scatter_point_dataset_->var_num);
        for (int i = 0; i < scatter_point_dataset_->var_num; ++i) {
            value_ranges[i].resize(2);
            value_ranges[i][0] = 1e20;
            value_ranges[i][1] = -1e20;
        }

		for (int i = 0; i < selection_index.size(); ++i) {
			parallel_dataset_->subset_names[i] = QString("Cluster %0").arg(i);
			parallel_dataset_->subset_colors[i] = QColor(cluster_color[3 * selection_index[i]], cluster_color[3 * selection_index[i] + 1], cluster_color[3 * selection_index[i] + 2]);
			for (int j = 0; j < scatter_point_dataset_->point_num; ++j)
				if (cluster_index[j] == selection_index[i]) {
					ParallelRecord* record = new ParallelRecord;
					record->values = scatter_point_dataset_->original_point_values[j];

                    for (int k = 0; k < record->values.size(); ++k) {
                        if (value_ranges[k][0] > record->values[k]) value_ranges[k][0] = record->values[k];
                        if (value_ranges[k][1] < record->values[k]) value_ranges[k][1] = record->values[k];
                    }

					parallel_dataset_->subset_records[i].push_back(record);
				}
		}

        for (int i = 0; i < scatter_point_dataset_->var_num; ++i)
            if (value_ranges[i][1] - value_ranges[i][0] == 0) {
                value_ranges[i][0] -= 0.5;
                value_ranges[i][1] += 0.5;
            }
        for (int i = 0; i < parallel_dataset_->subset_records.size(); ++i) {
            for (int j = 0; j < parallel_dataset_->subset_records[i].size(); ++j) {
                for (int k = 0; k < scatter_point_dataset_->var_num; ++k) {
                    parallel_dataset_->subset_records[i][j]->values[k] = (parallel_dataset_->subset_records[i][j]->values[k] - value_ranges[k][0]) / (value_ranges[k][1] - value_ranges[k][0]);
                }
            }
            
        }

		for (int i = 0; i < scatter_point_dataset_->original_value_ranges.size(); ++i) {
			parallel_dataset_->axis_names[i] = scatter_point_dataset_->var_names[i];
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(value_ranges[i][0]));
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(value_ranges[i][1]));
		}
		parallel_dataset_->CompleteInput();
		parallel_dataset_->UpdateGaussian();

		parallel_dataset_->is_updating = false;

		parallel_coordinate_->SetDataset(parallel_dataset_);
		Utility::GenerateAxisOrder(parallel_dataset_, var_axis_order);

        std::vector< int > focus_index;
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
        }

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
	if (transmap_data_ != NULL) {
		for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
			transmap_data_->cluster_nodes[i]->is_highlighted = false;
		transmap_data_->cluster_nodes.clear();
	}

	int cluster_num;
	std::vector< int > cluster_index;
	float dis_per_pixel = this->GetMainViewDisPerPixel();
	// dis_per_pixel * 100.0
	cluster_tree_->GetClusterResult(current_view_level_, transmap_data_->cluster_nodes);
	cluster_tree_->GetClusterResult(current_view_level_, cluster_num, cluster_index);
	
	std::vector< QColor > colors;
	for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
		colors.push_back(transmap_data_->cluster_nodes[i]->color);

	original_point_rendering_layer_->SetClusterIndex(cluster_num, cluster_index, colors);
	original_point_rendering_layer_->SetHighlightCluster(-1);

	transmap_data_->ClearData();
	transmap_data_->cluster_num = cluster_num;
	transmap_data_->var_num = scatter_point_dataset_->var_weights.size();
	transmap_data_->dataset = scatter_point_dataset_;
	transmap_data_->ProcessData();

	trans_map_->SetNodeRadius(dis_per_pixel * label_pixel_radius_);
	trans_map_->SetData(scatter_point_dataset_, transmap_data_);
	trans_map_->SetAxisOrder(var_axis_order);
    std::vector< QColor > point_colors;
    for (int i = 0; i < cluster_index.size(); ++i) {
        point_colors.push_back(colors[cluster_index[i]]);
    }
    trans_map_->UpdateDensityActor(point_colors);

	main_view_->update();
}


void ScatterPointGlyph::UpdatePointMap() {
	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}
}

void ScatterPointGlyph::UpdatePathMap() {

}

void ScatterPointGlyph::UpdateTreemap() {
	std::vector< int > selected_ids;
	trans_map_->GetSelectedClusterIds(selected_ids);

	//UncertaintyTree* un_tree = dynamic_cast<UncertaintyTree*>(cluster_tree_);
	cluster_tree_->SortTree(selected_ids);

	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	std::vector< bool > is_selected;
	is_selected.resize(transmap_data_->cluster_nodes.size(), false);

	std::vector< CNode* > selected_nodes;
	if (selected_ids.size() != 0) {
		for (int i = 0; i < selection_index.size(); ++i) {
			selected_nodes.push_back(transmap_data_->cluster_nodes[selection_index[i]]);
			is_selected[selection_index[i]] = true;
		}
		for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
			if (!is_selected[i]) selected_nodes.push_back(transmap_data_->cluster_nodes[i]);
		tree_map_view_->SetData(cluster_tree_, scatter_point_dataset_->var_num, selected_nodes, selected_ids.size(), var_axis_order, scatter_point_dataset_->var_names, scatter_point_dataset_->var_colors);
	}
	else {
		selected_nodes = transmap_data_->cluster_nodes;
		tree_map_view_->SetData(cluster_tree_, scatter_point_dataset_->var_num, selected_nodes, selected_nodes.size(), var_axis_order, scatter_point_dataset_->var_names, scatter_point_dataset_->var_colors);
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
		MultiLabelTree* un_tree = dynamic_cast< MultiLabelTree* >(cluster_tree_);

		int temp_cluster_index = trans_map_->GetSelectedClusterIndex();
		if (temp_cluster_index != -1) {
			un_tree->SplitCluster(transmap_data_->cluster_nodes[temp_cluster_index]->id());
			UpdateTransmap();
			UpdateTreemap();
			UpdateParallelCoordinate();
		}
	}
}

void ScatterPointGlyph::OnMergeClusterTriggered() {
	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* un_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);

		std::vector< int > merged_clusters;
		trans_map_->GetSelectedClusterIndex(merged_clusters);

		std::vector< int > cluster_ids;
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

	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}

	this->UpdateParallelCoordinate();

	this->main_view_->update();
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
	std::list< CNode* > node_seq = trans_map_->GetNodeSequence();

	PathRecord* record = new PathRecord;
	record->item_values.resize(node_seq.size());
	record->item_color.resize(node_seq.size());
	std::list< CNode* >::iterator node_iter = node_seq.begin();

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

		std::vector< float > values;
        for (int i = 0; i < scatter_point_dataset_->original_point_values.size(); ++i) 
            values.push_back(scatter_point_dataset_->original_point_values[i][var_index - 1]);
        original_point_rendering_layer_->SetPointValue(values);
        main_view_->update();
	} else {
		trans_map_->ShowVarTrend(-1);
        std::vector< float > values;
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
	label_pixel_radius_ = map_control_ui_.glyph_size_spinbox->value();

	this->UpdateTransmap();
}
