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

#include "layer_control_widget.h"
#include "rendering_layer_model.h"
#include "hier_para_widget.h"
#include "hier_solver.h"
#include "cluster_glyph_layer.h"
#include "continuity_extractor.h"
#include "point_rendering_layer.h"
#include "gestalt_processor2.h"
#include "scatter_point_dataset.h"
#include "grid_scatter_point_dataset.h"
#include "scatter_point_view.h"
#include "wrf_data_manager.h"
#include "map_rendering_layer.h"
#include "hierarchical_tree.h"
#include "immediate_tree.h"
#include "immediate_gestalt_tree.h"
#include "uncertainty_tree.h"
#include "parallel_coordinate.h"
#include "transmap.h"
#include "transmap_data.h"
#include "color_mapping_generator.h"
#include "field_rendering_layer.h"
#include "distance_matrix_layer.h"
#include "path_explore_widget.h"
#include "path_dataset.h"
#include "tour_path_generator.h"
#include "change_table_lens.h"
#include "tree_map_view.h"

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), dataset_(NULL), sys_mode_(UNCERTAINTY_MODE),
	original_point_rendering_layer_(NULL), cluster_point_rendering_layer_(NULL), un_rendering_layer_(NULL), 
	map_rendering_layer_(NULL), parallel_dataset_(NULL), dis_per_pixel_(0.0), min_pixel_radius_(5) {

	ui_.setupUi(this);

	data_manager_ = WrfDataManager::GetInstance();
	path_generator_ = new TourPathGenerator;
	pathset_ = new PathDataset;

	this->InitWidget();

	is_active_retrieval_on_ = false;

	cluster_tree_vec_.resize(5, NULL);
}

ScatterPointGlyph::~ScatterPointGlyph() {

}

void ScatterPointGlyph::InitWidget() {
	QList< int > temp_sizes;
	temp_sizes.push_back(200);
	temp_sizes.push_back(100);

	main_view_ = new ScatterPointView;
	main_view_->setFocusPolicy(Qt::StrongFocus);
	parallel_coordinate_ = new ParallelCoordinate;
	QSplitter* transmap_widget_splitter = new QSplitter(Qt::Vertical);
	transmap_widget_splitter->addWidget(main_view_);
	transmap_widget_splitter->addWidget(parallel_coordinate_);
	transmap_widget_splitter->setStretchFactor(0, 2);
	transmap_widget_splitter->setStretchFactor(1, 1);

	this->setDockOptions(QMainWindow::AllowNestedDocks);

	tree_map_view_ = new TreeMapView;
	tree_map_view_->setFixedWidth(700);
	tree_map_panel_ = new QDockWidget(QString("Tree Map Panel"), this);
	tree_map_panel_->setWidget(tree_map_view_);
	this->addDockWidget(Qt::RightDockWidgetArea, tree_map_panel_);
	ui_.menuView->addAction(tree_map_panel_->toggleViewAction());
	tree_map_panel_->setVisible(false);

	path_explore_view_ = new PathExploreWidget;
	path_explore_view_->setFixedWidth(700);
	path_explore_panel_ = new QDockWidget(QString("Path Explore Panel"), this);
	path_explore_panel_->setWidget(path_explore_view_);
	this->tabifyDockWidget(tree_map_panel_, path_explore_panel_);
	ui_.menuView->addAction(path_explore_panel_->toggleViewAction());
	path_explore_panel_->setVisible(false);

	QHBoxLayout* main_layout = new QHBoxLayout;
	main_layout->addWidget(transmap_widget_splitter);
	ui_.centralWidget->setLayout(main_layout);

	main_renderer_ = vtkRenderer::New();
	main_view_->GetRenderWindow()->AddRenderer(main_renderer_);
	main_renderer_->SetViewport(0.0, 0.0, 1.0, 1.0);
	main_renderer_->SetBackground(1.0, 1.0, 1.0);
	connect(main_view_, SIGNAL(ViewUpdated()), this, SLOT(OnMainViewUpdated()));
	connect(main_view_, SIGNAL(GlyphSelected(int, int)), this, SLOT(OnGlyphSelected(int, int)));
	connect(main_view_, SIGNAL(LeftButtonUp()), this, SLOT(OnMainviewLeftButtonUp()));
	connect(main_view_, SIGNAL(RightButtonDown()), this, SLOT(OnMainViewRightButtonDown()));
	connect(main_view_, SIGNAL(MouseDrag(int, int)), this, SLOT(OnMouseDragmove(int, int)));

	/*dis_matrix_renderer_ = vtkRenderer::New();
	main_view_->GetRenderWindow()->AddRenderer(dis_matrix_renderer_);
	dis_matrix_renderer_->SetViewport(0.7, 0.5, 1.0, 1.0);
	double eyepos[] = { 0.5, 0.5, 2.0 };
	double viewup[] = { 0.0, 1.0, 0.0 };
	dis_matrix_renderer_->GetActiveCamera()->SetPosition(eyepos);
	dis_matrix_renderer_->GetActiveCamera()->SetFocalPoint(0.5, 0.5, 0.0);
	dis_matrix_renderer_->GetActiveCamera()->Modified();
	dis_matrix_renderer_->SetBackground(1.0, 1.0, 1.0);

	other_renderer_ = vtkRenderer::New();
	main_view_->GetRenderWindow()->AddRenderer(other_renderer_);
	other_renderer_->SetViewport(0.7, 0.0, 1.0, 0.5);
	other_renderer_->SetBackground(1.0, 1.0, 1.0);*/

	layer_control_widget_ = new LayerControlWidget;
	rendering_layer_model_ = new RenderingLayerModel;
	layer_control_widget_->SetRenderingLayerModel(rendering_layer_model_);
	layer_control_panel_ = new QDockWidget(QString("Layer Control Panel"), this);
	layer_control_panel_->setWidget(layer_control_widget_);
	this->tabifyDockWidget(tree_map_panel_, layer_control_panel_);
	ui_.menuView->addAction(layer_control_panel_->toggleViewAction());
	layer_control_panel_->setVisible(false);
	connect(rendering_layer_model_, SIGNAL(LayerPropertyChanged()), main_view_, SLOT(update()));

	trans_map_ = new TransMap(main_view_);
	trans_map_->SetInteractor(main_view_->GetInteractor());
	trans_map_->SetDefaultRenderer(this->main_renderer_);
	transmap_data_ = new TransMapData;
	rendering_layer_model_->AddLayer(QString("Transfer Map"), trans_map_, false);

	hier_para_widget_ = new HierParaWidget;
	hier_para_panel_ = new QDockWidget(QString("Hierarchical Clustering Panel"), this);
	hier_para_panel_->setWidget(hier_para_widget_);
	hier_para_panel_->setVisible(false);
	this->tabifyDockWidget(tree_map_panel_, hier_para_panel_);
	ui_.menuView->addAction(hier_para_panel_->toggleViewAction());

	sys_mode_action_group_ = new QActionGroup(this);
	sys_mode_action_group_->addAction(ui_.action_perception_driven);
	sys_mode_action_group_->addAction(ui_.action_hierarchical_clustering);
	sys_mode_action_group_->addAction(ui_.actionImmediate_Gestalts);
	sys_mode_action_group_->setExclusive(true);

	connect(ui_.actionVTK_Unstructured_Grid_Data, SIGNAL(triggered()), this, SLOT(OnActionOpenScatterFileTriggered()));
	connect(ui_.actionRaw_Grid_Data, SIGNAL(triggered()), this, SLOT(OnActionOpenRawGridFileTriggered()));
	connect(ui_.actionExec, SIGNAL(triggered()), this, SLOT(OnExecClusteringTriggered()));
	connect(ui_.actionBrush_Cluster, SIGNAL(toggled(bool)), this, SLOT(OnBrushSelectionTriggered(bool)));
	connect(ui_.actionSplit, SIGNAL(triggered()), this, SLOT(OnSplitClusterTriggered()));
	connect(ui_.actionMerge, SIGNAL(triggered()), this, SLOT(OnMergeClusterTriggered()));
	connect(ui_.actionAdd_Path_Sequence, SIGNAL(triggered()), this, SLOT(OnAddPathSequenceTriggered()));
}

void ScatterPointGlyph::OnActionOpenVtkFileTriggered() {
	if (dataset_ == NULL) dataset_ = new ScatterPointDataset;

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

			dataset_->original_point_pos.push_back(pos);
			dataset_->original_point_values.push_back(value);
		}
	}

	dataset_->DirectConstruct();
	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionOpenRawGridFileTriggered() {
	GridScatterPointDataset* grid_dataset = new GridScatterPointDataset;

	data_manager_->LoadEnsembleData(WRF_ACCUMULATED_PRECIPITATION, std::string("E:/Data/ens_l/apcp_sfc_latlon_all_20150401_20150430_liaoLeeVC4.nc"));
	data_manager_->LoadEnsembleData(WRF_MSLP, std::string("E:/Data/ens_l/pres_msl_latlon_all_20150401_20150430_liaoawNRCR.nc"));
	data_manager_->LoadEnsembleData(WRF_T2M, std::string("E:/Data/ens_l/tmp_2m_latlon_all_20150401_20150430_liaoSwbPtU.nc"));
	data_manager_->LoadEnsembleData(WRF_PRECIPITABLE_WATER, std::string("E:/Data/ens_l/pwat_eatm_latlon_all_20150401_20150430_liao5WubTq.nc"));
	std::vector< std::vector< WrfGridValueMap* > > maps;
	maps.resize(4);
	QDateTime temp_datetime = QDateTime(QDate(2015, 4, 5), QTime(0, 0));
	data_manager_->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_ACCUMULATED_PRECIPITATION, 24, maps[0]);
	data_manager_->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_MSLP, 24, maps[1]);
	data_manager_->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_T2M, 24, maps[2]);
	data_manager_->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_PRECIPITABLE_WATER, 24, maps[3]);
	MapRange map_range;
	data_manager_->GetModelMapRange(WRF_NCEP_ENSEMBLES, map_range);

	int point_num = map_range.x_grid_number * map_range.y_grid_number;
	grid_dataset->original_point_values.resize(point_num);
	grid_dataset->original_point_pos.resize(point_num);
	for (int i = 0; i < point_num; ++i) grid_dataset->original_point_pos[i].resize(2);

	for (int i = 0; i < point_num; ++i){
		int w = i % map_range.x_grid_number;
		int h = i / map_range.x_grid_number;
		grid_dataset->original_point_pos[i][0] = map_range.start_x + map_range.x_grid_space * w;
		grid_dataset->original_point_pos[i][1] = map_range.start_y + map_range.y_grid_space * h;

		grid_dataset->original_point_values[i].push_back(maps[0][0]->values[i]);
		grid_dataset->original_point_values[i].push_back(maps[1][0]->values[i]);
		grid_dataset->original_point_values[i].push_back(maps[2][0]->values[i]);
		grid_dataset->original_point_values[i].push_back(maps[3][0]->values[i]);
	}
	grid_dataset->w = map_range.x_grid_number;
	grid_dataset->h = map_range.y_grid_number;
	grid_dataset->DirectConstruct();

	if (dataset_ != NULL) delete dataset_;
	dataset_ = grid_dataset;

	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionOpenScatterFileTriggered() {
	if (dataset_ == NULL) dataset_ = new ScatterPointDataset;

	std::ifstream input_file("E:/Projects/DataProcessing/MdsConverter/mds_result.txt");
	int record_num, value_num;
	input_file >> record_num >> value_num;
	dataset_->original_point_pos.resize(record_num);
	dataset_->original_point_values.resize(record_num);
	for (int i = 0; i < record_num; ++i) {
		dataset_->original_point_pos[i].resize(2);
		dataset_->original_point_values[i].resize(7);

		float temp;
		input_file >> dataset_->original_point_pos[i][0] >> dataset_->original_point_pos[i][1];
		for (int j = 0; j < value_num; ++j)
			if (j < 7)
				input_file >> dataset_->original_point_values[i][j];
			else
				input_file >> temp;
	}

	/*std::ifstream input_file("E:/Projects/DataProcessing/MdsConverter/dingpao_1.txt");
	int record_num, value_num;
	input_file >> record_num >> value_num;
	dataset_->original_point_pos.resize(record_num);
	dataset_->original_point_values.resize(record_num);
	for (int i = 0; i < record_num; ++i) {
		dataset_->original_point_pos[i].resize(2);
		dataset_->original_point_values[i].resize(value_num);

		float temp;
		input_file >> dataset_->original_point_pos[i][0] >> dataset_->original_point_pos[i][1];
		for (int j = 0; j < value_num; ++j)
			if (j < 7)
				input_file >> dataset_->original_point_values[i][j];
			else
				input_file >> temp;
	}*/

	dataset_->DirectConstruct();
	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionCloseTriggered() {

}

void ScatterPointGlyph::OnActionExitTriggered() {
	
}

void ScatterPointGlyph::OnSysmodeChanged() {
	if (ui_.action_perception_driven->isChecked()) {
		sys_mode_ = PERCEPTION_MODE;
	} else if (ui_.action_hierarchical_clustering->isChecked()) {
		sys_mode_ = HIER_MODE;
	} else if (ui_.actionImmediate_Gestalts->isChecked()) {
		sys_mode_ = IMMEDIATE_PERCEPTION_MODE;
	}
}

void ScatterPointGlyph::OnExecClusteringTriggered() {
	if (dataset_ == NULL) {
		QMessageBox::information(this, tr("Warning"), tr("Please load data first."));
		return;
	}

	this->dis_per_pixel_ = this->GetMainViewDisPerPixel();
	
	switch (sys_mode_)
	{
	case ScatterPointGlyph::HIER_MODE:
	{
		if (cluster_tree_vec_[HIER_MODE] == NULL) {
			cluster_tree_vec_[HIER_MODE] = new HierarchicalTree(dataset_);

			connect(cluster_tree_vec_[HIER_MODE], SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}
	}
		break;
	case ScatterPointGlyph::PERCEPTION_MODE:
	{
		if (cluster_tree_vec_[PERCEPTION_MODE] == NULL) {
			cluster_tree_vec_[PERCEPTION_MODE] = new HierarchicalTree(dataset_);

			connect(cluster_tree_vec_[PERCEPTION_MODE], SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}

		HierarchicalTree* hier_tree = dynamic_cast< HierarchicalTree* >(cluster_tree_vec_[PERCEPTION_MODE]);
		if (hier_tree != NULL) hier_tree->start();
	}
		break;
	case ScatterPointGlyph::IMMEDIATE_PERCEPTION_MODE:
	{
		if (cluster_tree_vec_[IMMEDIATE_PERCEPTION_MODE] == NULL) {
			cluster_tree_vec_[IMMEDIATE_PERCEPTION_MODE] = new ImmediateTree(dataset_);

			connect(cluster_tree_vec_[IMMEDIATE_PERCEPTION_MODE], SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}

		ImmediateTree* immediate_tree = dynamic_cast< ImmediateTree* >(cluster_tree_vec_[IMMEDIATE_PERCEPTION_MODE]);
		if (immediate_tree != NULL) {
			float left, right, bottom, top;
			this->GetSceneRange(left, right, bottom, top);
			dataset_->Sample(left, right, bottom, top);
			immediate_tree->SetSampleSize((int)(this->main_view_->width() * this->main_view_->height() / (80 * 80)));
			immediate_tree->SetRadiusThreshold(100.0 * this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]));
			immediate_tree->start();
		}
	}
		break;
	case ScatterPointGlyph::IMMEDIATE_GESTALT_MODE:
	{
		if (cluster_tree_vec_[IMMEDIATE_GESTALT_MODE] == NULL) {
			cluster_tree_vec_[IMMEDIATE_GESTALT_MODE] = new ImmediateGestaltTree(dataset_);

			connect(cluster_tree_vec_[IMMEDIATE_GESTALT_MODE], SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}

		ImmediateGestaltTree* immediate_tree = dynamic_cast< ImmediateGestaltTree* >(cluster_tree_vec_[IMMEDIATE_GESTALT_MODE]);
		if (immediate_tree != NULL) {
			float left, right, bottom, top;
			this->GetSceneRange(left, right, bottom, top);
			dataset_->Sample(left, right, bottom, top);
			immediate_tree->SetSampleSize((int)(this->main_view_->width() * this->main_view_->height() / (160 * 160)));
			immediate_tree->SetRadiusThreshold(100.0 * this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]));
			immediate_tree->start();
		}
	}
		break;
	case ScatterPointGlyph::UNCERTAINTY_MODE:
	{
		if (cluster_tree_vec_[UNCERTAINTY_MODE] == NULL) {
			cluster_tree_vec_[UNCERTAINTY_MODE] = new UncertaintyTree(dataset_);

			connect(cluster_tree_vec_[UNCERTAINTY_MODE], SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}

		UncertaintyTree* un_tree = dynamic_cast< UncertaintyTree* >(cluster_tree_vec_[UNCERTAINTY_MODE]);
		if (un_tree != NULL) {
			float left, right, bottom, top;
			this->GetSceneRange(left, right, bottom, top);
			dataset_->Sample(left, right, bottom, top);
			un_tree->SetRadiusThreshold(200.0 * this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]));
			un_tree->start();
		}
	}
		break;
	default:
		break;
	}
}

void ScatterPointGlyph::AddPointData2View() {
	if (original_point_rendering_layer_ == NULL) {
		original_point_rendering_layer_ = new PointRenderingLayer;
		original_point_rendering_layer_->SetInteractor(main_view_->GetInteractor());
		original_point_rendering_layer_->SetDefaultRenderer(this->main_renderer_);
		rendering_layer_model_->AddLayer(QString("Original Point Layer"), original_point_rendering_layer_);
	}
	original_point_rendering_layer_->SetData(dataset_);

	main_renderer_->ResetCamera();

	this->UpdateParallelCoordinate();

	main_view_->update();
}

void ScatterPointGlyph::UpdateParallelCoordinate() {
	if (parallel_dataset_ == NULL) {
		parallel_dataset_ = new ParallelDataset;
		parallel_dataset_->subset_names.push_back(QString("MDS Data"));
		parallel_dataset_->subset_records.resize(1);
		parallel_dataset_->axis_anchors.resize(dataset_->original_point_values[0].size());
		for (int i = 0; i < dataset_->original_point_values[0].size(); ++i) {
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(dataset_->original_value_ranges[i][0]));
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(dataset_->original_value_ranges[i][1]));
			parallel_dataset_->axis_names.push_back(QString("v"));
		}

		for (int i = 0; i < dataset_->point_values.size(); ++i) {
			ParallelRecord* record = new ParallelRecord;
			record->values = dataset_->point_values[i];
			parallel_dataset_->subset_records[0].push_back(record);
		}
		parallel_dataset_->CompleteInput();
		parallel_coordinate_->SetDataset(parallel_dataset_);

		QString str = QString("%0 points.").arg(dataset_->original_point_pos.size());
		ui_.statusBar->showMessage(str);

	} else {
		std::vector< int > selection_index;
		trans_map_->GetSelectedClusterIndex(selection_index);
		this->GenerateParallelDataset(parallel_dataset_, selection_index);
		parallel_coordinate_->update();

		int selected_count = 0;
		for (int i = 0; i < parallel_dataset_->subset_records.size(); ++i)
			selected_count += parallel_dataset_->subset_records[i].size();
		QString str = QString("%0 out of %1 (%2%) points are selected.").arg(selected_count).arg(dataset_->original_point_pos.size()).arg((int)selected_count * 100 / dataset_->original_point_pos.size());
		ui_.statusBar->showMessage(str);
	}
}

void ScatterPointGlyph::GenerateParallelDataset(ParallelDataset* pdata, std::vector< int >& cluster_ids) {
	pdata->is_updating = true;
	pdata->ClearData();

	std::vector< int > cluster_color = original_point_rendering_layer_->GetClusterColor();
	parallel_dataset_->subset_names.resize(cluster_ids.size());
	parallel_dataset_->subset_colors.resize(cluster_ids.size());
	parallel_dataset_->subset_records.resize(cluster_ids.size());
	parallel_dataset_->axis_names.resize(dataset_->original_value_ranges.size());
	parallel_dataset_->axis_anchors.resize(dataset_->original_value_ranges.size());

	for (int i = 0; i < cluster_ids.size(); ++i) {
		parallel_dataset_->subset_names[i] = QString("Cluster %0").arg(cluster_ids[i]);
		parallel_dataset_->subset_colors[i] = QColor(cluster_color[3 * cluster_ids[i]], cluster_color[3 * cluster_ids[i] + 1], cluster_color[3 * cluster_ids[i] + 2]);
		for (int j = 0; j < dataset_->sample_index.size(); ++j)
			if (cluster_index[dataset_->sample_index[j]] == cluster_ids[i]) {
				ParallelRecord* record = new ParallelRecord;
				record->values = dataset_->point_values[j];
				parallel_dataset_->subset_records[i].push_back(record);
			}
	}
	for (int i = 0; i < dataset_->original_value_ranges.size(); ++i) {
		parallel_dataset_->axis_names[i] = "v";
		parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(dataset_->original_value_ranges[i][0]));
		parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(dataset_->original_value_ranges[i][1]));
	}
	parallel_dataset_->CompleteInput();
	parallel_dataset_->UpdateGaussian();

	pdata->is_updating = false;
}

void ScatterPointGlyph::OnClusterFinished() {
	if (sys_mode_ == PERCEPTION_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		TreeCommon* tree = cluster_tree_vec_[PERCEPTION_MODE];
		tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluster_num, cluster_index);
	}

	if (sys_mode_ == IMMEDIATE_PERCEPTION_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		TreeCommon* tree = cluster_tree_vec_[IMMEDIATE_PERCEPTION_MODE];
		tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluster_num, cluster_index);
	}

	if (sys_mode_ == IMMEDIATE_GESTALT_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		TreeCommon* tree = cluster_tree_vec_[IMMEDIATE_GESTALT_MODE];
		tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluster_num, cluster_index);
	}

	if (sys_mode_ == UNCERTAINTY_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		UpdateTransmap();
		UpdateParallelCoordinate();
		UpdateTreemap();
		UpdatePathMap();
	}
}

void ScatterPointGlyph::OnMainViewUpdated() {
	this->UpdateTransmap();
	this->UpdateTreemap();
}

void ScatterPointGlyph::UpdateTransmap() {
	if (transmap_data_ != NULL) {
		for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
			transmap_data_->cluster_nodes[i]->is_highlighted = false;
	}

	UncertaintyTree* un_tree = dynamic_cast<UncertaintyTree*>(cluster_tree_vec_[UNCERTAINTY_MODE]);
	this->dis_per_pixel_ = this->GetMainViewDisPerPixel();
	if (is_active_retrieval_on_) {
		un_tree->GetActiveClusterResult(transmap_data_->cluster_nodes);
		un_tree->GetActiveClusterResult(cluster_num, cluster_index);
	} else {
		un_tree->GetClusterResult(this->dis_per_pixel_ * 100.0, transmap_data_->cluster_nodes);
		un_tree->GetClusterResult(this->dis_per_pixel_ * 100.0, cluster_num, cluster_index);
	}
	
	std::vector< QColor > colors;
	for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
		colors.push_back(transmap_data_->cluster_nodes[i]->color);

	original_point_rendering_layer_->SetClusterIndex(cluster_num, cluster_index, colors);
	original_point_rendering_layer_->SetHighlightCluster(-1);

	transmap_data_->ClearData();
	transmap_data_->cluster_num = cluster_num;
	transmap_data_->var_num = dataset_->weights.size();
	transmap_data_->dataset = dataset_;
	transmap_data_->ProcessData();

	trans_map_->SetNodeRadius(dis_per_pixel_ * 30);
	trans_map_->SetOriginalData(dataset_);
	trans_map_->SetData(transmap_data_);

	main_view_->update();
}


void ScatterPointGlyph::UpdatePathMap() {

}

void ScatterPointGlyph::UpdateTreemap() {
	UncertaintyTree* un_tree = dynamic_cast<UncertaintyTree*>(cluster_tree_vec_[UNCERTAINTY_MODE]);
	tree_map_view_->SetData(un_tree->root());
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

void ScatterPointGlyph::OnBrushSelectionTriggered(bool checked) {
	if (checked)
		trans_map_->SetBrushSelectionOn();
	else {
		trans_map_->SetBrushSelectionOff();
		this->main_view_->update();
	}
}

void ScatterPointGlyph::OnGlyphSelected(int x, int y) {
	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}

	UpdateParallelCoordinate();
	UpdateTreemap();
}

void ScatterPointGlyph::OnMainviewLeftButtonUp() {
	if (trans_map_ != NULL) {
		trans_map_->SetMouseReleased();
		if (ui_.actionBrush_Cluster->isChecked()) {
			std::vector< int > selection_index;
			trans_map_->GetSelectedClusterIndex(selection_index);
			if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
				original_point_rendering_layer_->SetHighlightClusters(selection_index);
			}

			UpdateParallelCoordinate();
			UpdateTreemap();

			this->main_view_->update();
		}
	}
}

void ScatterPointGlyph::OnMainViewRightButtonDown() {
	if (trans_map_->IsMapUpdateNeeded()) {
		UpdateTransmap();
		tree_map_view_->scene()->update();
	}
}

void ScatterPointGlyph::OnMouseDragmove(int x, int y) {
	if (trans_map_ != NULL) trans_map_->SetMouseDragmove(x, y);
}

void ScatterPointGlyph::OnSplitClusterTriggered() {
	if (sys_mode_ == UNCERTAINTY_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		UncertaintyTree* un_tree = dynamic_cast< UncertaintyTree* >(cluster_tree_vec_[UNCERTAINTY_MODE]);

		int temp_cluster_index = trans_map_->GetSelectedClusterIndex();
		if (temp_cluster_index != -1) {
			un_tree->SplitCluster(transmap_data_->cluster_nodes[temp_cluster_index]->id);
			UpdateTransmap();
			UpdateTreemap();
			UpdateParallelCoordinate();
		}
	}
}

void ScatterPointGlyph::OnMergeClusterTriggered() {
	if (sys_mode_ == UNCERTAINTY_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		UncertaintyTree* un_tree = dynamic_cast<UncertaintyTree*>(cluster_tree_vec_[UNCERTAINTY_MODE]);

		std::vector< int > merged_clusters;
		trans_map_->GetSelectedClusterIndex(merged_clusters);

		std::vector< int > cluster_ids;
		for (int i = 0; i < merged_clusters.size(); ++i)
			cluster_ids.push_back(transmap_data_->cluster_nodes[merged_clusters[i]]->id);
		un_tree->MergeClusters(cluster_ids);

		UpdateTransmap();
		UpdateTreemap();
		UpdateParallelCoordinate();
	}
}

void ScatterPointGlyph::OnTreemapNodeSelected(int node_id) {
	trans_map_->SetNodeSelected(node_id);

	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}

	this->GenerateParallelDataset(parallel_dataset_, selection_index);

	parallel_coordinate_->update();

	this->main_view_->update();
}

void ScatterPointGlyph::OnAddPathSequenceTriggered() {
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

	pathset_->path_records.push_back(record);
	path_explore_view_->SetData(pathset_);
}