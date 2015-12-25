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

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), dataset_(NULL), sys_mode_(UNCERTAINTY_MODE),
	original_point_rendering_layer_(NULL), cluster_point_rendering_layer_(NULL), un_rendering_layer_(NULL), 
	map_rendering_layer_(NULL), parallel_dataset_(NULL), dis_per_pixel_(0.0), min_pixel_radius_(5) {

	ui_.setupUi(this);

	data_manager_ = WrfDataManager::GetInstance();

	this->InitWidget();

	cluster_tree_vec_.resize(5, NULL);
}

ScatterPointGlyph::~ScatterPointGlyph() {

}

void ScatterPointGlyph::InitWidget() {
	main_view_ = new ScatterPointView;
	main_view_->setFocusPolicy(Qt::StrongFocus);
	parallel_coordinate_ = new ParallelCoordinate;
	QVBoxLayout* transmap_widget_layout = new QVBoxLayout;
	transmap_widget_layout->addWidget(main_view_);
	transmap_widget_layout->addWidget(parallel_coordinate_);

	path_explore_view_ = new PathExploreWidget;
	QVBoxLayout* path_explore_layout = new QVBoxLayout;
	path_explore_layout->addWidget(path_explore_view_);
	QWidget* temp_widget = new QWidget;
	path_explore_layout->addWidget(temp_widget);

	QHBoxLayout* main_layout = new QHBoxLayout;
	main_layout->addLayout(transmap_widget_layout);
	main_layout->addLayout(path_explore_layout);
	ui_.centralWidget->setLayout(main_layout);

	main_renderer_ = vtkRenderer::New();
	main_view_->GetRenderWindow()->AddRenderer(main_renderer_);
	main_renderer_->SetViewport(0.0, 0.0, 0.7, 1.0);
	main_renderer_->SetBackground(1.0, 1.0, 1.0);
	connect(main_view_, SIGNAL(ViewUpdated()), this, SLOT(UpdateClusterView()));
	connect(main_view_, SIGNAL(GlyphSelected(int, int)), this, SLOT(OnGlyphSelected(int, int)));

	dis_matrix_renderer_ = vtkRenderer::New();
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
	//other_renderer_->SetBackground(1.0, 1.0, 1.0);

	layer_control_widget_ = new LayerControlWidget;
	rendering_layer_model_ = new RenderingLayerModel;
	layer_control_widget_->SetRenderingLayerModel(rendering_layer_model_);
	layer_control_panel_ = new QDockWidget(QString("Layer Control Panel"), this);
	layer_control_panel_->setWidget(layer_control_widget_);
	this->addDockWidget(Qt::RightDockWidgetArea, layer_control_panel_);
	ui_.menuView->addAction(layer_control_panel_->toggleViewAction());
	connect(rendering_layer_model_, SIGNAL(LayerPropertyChanged()), main_view_, SLOT(update()));

	trans_map_ = TransMap::New();
	trans_map_->SetInteractor(main_view_->GetInteractor());
	trans_map_->SetDefaultRenderer(this->main_renderer_);
	rendering_layer_model_->AddLayer(QString("Transfer Map"), trans_map_, false);

	hier_para_widget_ = new HierParaWidget;
	hier_para_panel_ = new QDockWidget(QString("Hierarchical Clustering Panel"), this);
	hier_para_panel_->setWidget(hier_para_widget_);
	hier_para_panel_->setVisible(false);
	this->addDockWidget(Qt::RightDockWidgetArea, hier_para_panel_);
	ui_.menuView->addAction(hier_para_panel_->toggleViewAction());

	sys_mode_action_group_ = new QActionGroup(this);
	sys_mode_action_group_->addAction(ui_.action_perception_driven);
	sys_mode_action_group_->addAction(ui_.action_hierarchical_clustering);
	sys_mode_action_group_->addAction(ui_.actionImmediate_Gestalts);
	sys_mode_action_group_->setExclusive(true);

	connect(ui_.actionVTK_Unstructured_Grid_Data, SIGNAL(triggered()), this, SLOT(OnActionOpenScatterFileTriggered()));
	connect(ui_.actionRaw_Grid_Data, SIGNAL(triggered()), this, SLOT(OnActionOpenRawGridFileTriggered()));
	connect(ui_.actionExec, SIGNAL(triggered()), this, SLOT(UpdateClusterView()));
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
	if (dataset_ == NULL) dataset_ = new ScatterPointDataset;

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
	dataset_->original_point_values.resize(point_num);
	dataset_->original_point_pos.resize(point_num);
	for (int i = 0; i < point_num; ++i) dataset_->original_point_pos[i].resize(2);

	for (int i = 0; i < point_num; ++i){
		int w = i % map_range.x_grid_number;
		int h = i / map_range.x_grid_number;
		dataset_->original_point_pos[i][0] = map_range.start_x + map_range.x_grid_space * w;
		dataset_->original_point_pos[i][1] = map_range.start_y + map_range.y_grid_space * h;

		dataset_->original_point_values[i].push_back(maps[0][0]->values[i]);
		dataset_->original_point_values[i].push_back(maps[1][0]->values[i]);
		dataset_->original_point_values[i].push_back(maps[2][0]->values[i]);
		dataset_->original_point_values[i].push_back(maps[3][0]->values[i]);
	}
	dataset_->is_structured_data = true;
	dataset_->w = map_range.x_grid_number;
	dataset_->h = map_range.y_grid_number;

	dataset_->DirectConstruct();
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
		dataset_->original_point_values[i].resize(value_num);

		input_file >> dataset_->original_point_pos[i][0] >> dataset_->original_point_pos[i][1];
		for (int j = 0; j < value_num; ++j)
			input_file >> dataset_->original_point_values[i][j];
	}

	dataset_->DirectConstruct();
	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionExtractDataTriggered() {
	
}

void ScatterPointGlyph::OnActionCloseTriggered() {

}

void ScatterPointGlyph::OnActionExitTriggered() {
	
}

void ScatterPointGlyph::OnGlyphSelected(int x, int y) {
	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}

	this->GenerateParallelDataset(parallel_dataset_, selection_index);

	parallel_coordinate_->update();
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

void ScatterPointGlyph::UpdateClusterView() {
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
			un_tree->SetSampleSize((int)(this->main_view_->width() * this->main_view_->height() / (160 * 160)));
			un_tree->SetRadiusThreshold(100.0 * this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]));
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

	// For test
	for (int i = 0; i < dataset_->weights.size(); ++i) {
		FieldRenderingLayer* field_layer = new FieldRenderingLayer;
		field_layer->SetFieldData(dataset_, i);
		std::vector< ColorMappingGenerator::ColorIndex > color_index;
		bool b;
		ColorMappingGenerator::GetInstance()->GetColorIndex(RMS_MAPPING, color_index, b);
		std::vector< float > values;
		std::vector< double > rgb;
		for (int i = 0; i < color_index.size(); ++i) {
			values.push_back(color_index[i].value_index);
			rgb.push_back(color_index[i].color.redF());
			rgb.push_back(color_index[i].color.greenF());
			rgb.push_back(color_index[i].color.blueF());
			rgb.push_back(1.0);
		}
		for (int j = 0; j < values.size(); ++j) values[j] /= values[values.size() - 1];
		field_layer->SetColorScalar(values, rgb);
		field_layer->SetInteractor(this->main_view_->GetInteractor());
		field_layer->SetDefaultRenderer(this->main_renderer_);
		rendering_layer_model_->AddLayer(QString("Field Layer"), field_layer, false);
	}

	DistanceMatrixLayer* matrix = new DistanceMatrixLayer;
	std::vector< std::string > names;
	std::vector< std::vector< float > > distance;
	for (int i = 0; i < 5; ++i) names.push_back(std::string("test"));
	distance.resize(5);
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) distance[i].push_back((float)rand() / RAND_MAX);
	}
	matrix->SetData(names, distance);
	matrix->SetInteractor(this->main_view_->GetInteractor());
	matrix->SetDefaultRenderer(this->dis_matrix_renderer_);
	matrix->SetEnabled(true);

	this->UpdateParallelCoordinate();

	// for test
	PathDataset* pathset = new PathDataset;
	for (int i = 0; i < 5; ++i) {
		PathRecord* record = new PathRecord;
		record->item_values.resize(4);
		for (int j = 0; j < 4; ++j) record->item_values[j].resize(18, 1);
		record->change_values.resize(3);
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 18; ++k)
				record->change_values[j].push_back((float)rand() / RAND_MAX - 0.5);
		}
		pathset->path_records.push_back(record);
		for (int k = 0; k < 18; ++k) record->var_names.push_back("Test");
	}
	path_explore_view_->SetData(pathset);

	main_view_->update();
}

void ScatterPointGlyph::UpdateParallelCoordinate() {
	if (parallel_dataset_ == NULL) {
		parallel_dataset_ = new ParallelDataset;
		parallel_dataset_->subset_names.push_back(QString("MDS Data"));
		parallel_dataset_->subset_records.resize(1);
		parallel_dataset_->axis_anchors.resize(dataset_->original_point_values[0].size());
		for (int i = 0; i < dataset_->original_point_values[0].size(); ++i) {
			parallel_dataset_->axis_anchors[i].push_back(QString("0"));
			parallel_dataset_->axis_anchors[i].push_back(QString("1"));
			parallel_dataset_->axis_names.push_back(QString("v"));
		}

		for (int i = 0; i < dataset_->point_values.size(); ++i) {
			ParallelRecord* record = new ParallelRecord;
			record->values = dataset_->point_values[i];
			parallel_dataset_->subset_records[0].push_back(record);
		}
		parallel_dataset_->CompleteInput();
		parallel_coordinate_->SetDataset(parallel_dataset_);
	}
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
		UncertaintyTree* un_tree = dynamic_cast< UncertaintyTree* >(cluster_tree_vec_[UNCERTAINTY_MODE]);
		un_tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluster_num, cluster_index);

		original_point_rendering_layer_->SetClusterIndex(cluster_num, cluster_index);

		TransMapData* data = new TransMapData;
		data->cluster_num = cluster_num;
		data->var_num = dataset_->weights.size();
		data->cluster_index = cluster_index;
		// TODO: update the cluster color
		data->cluster_reprentative_color = original_point_rendering_layer_->GetClusterColor();
		// TODO: update the variable color
		data->var_repsentative_color.resize(data->var_num * 3, 0);
		data->dataset = dataset_;
		data->ProcessData();

		trans_map_->SetNodeRadius(dis_per_pixel_ * 20);
		trans_map_->SetData(data);
		trans_map_->SetEnabled(true);

		layer_control_widget_->UpdateWidget();
	}

	main_view_->update();
}

void ScatterPointGlyph::GenerateParallelDataset(ParallelDataset* pdata, std::vector< int >& cluster_ids) {
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
		for (int j = 0; j < cluster_index.size(); ++j)
			if (cluster_index[j] == cluster_ids[i]) {
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