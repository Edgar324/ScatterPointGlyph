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

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), dataset_(NULL), sys_mode_(UNCERTAINTY_MODE),
	cluster_glyph_layer_(NULL), original_point_rendering_layer_(NULL), cluster_point_rendering_layer_(NULL), un_rendering_layer_(NULL), 
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
	QVBoxLayout* central_widget_layout = new QVBoxLayout;
	central_widget_layout->addWidget(main_view_);
	central_widget_layout->addWidget(parallel_coordinate_);
	ui_.centralWidget->setLayout(central_widget_layout);

	trans_map_ = TransMap::New();
	trans_map_->SetInteractor(main_view_->GetInteractor());

	main_renderer_ = vtkRenderer::New();
	main_renderer_->RemoveAllViewProps();
	main_renderer_->SetBackground(1.0, 1.0, 1.0);
	main_view_->GetRenderWindow()->AddRenderer(main_renderer_);
	connect(main_view_, SIGNAL(ViewUpdated()), this, SLOT(UpdateClusterView()));
	connect(main_view_, SIGNAL(GlyphSelected(int, int)), this, SLOT(OnGlyphSelected(int, int)));

	layer_control_widget_ = new LayerControlWidget;
	rendering_layer_model_ = new RenderingLayerModel;
	layer_control_widget_->SetRenderingLayerModel(rendering_layer_model_);
	layer_control_panel_ = new QDockWidget(QString("Layer Control Panel"), this);
	layer_control_panel_->setWidget(layer_control_widget_);
	this->addDockWidget(Qt::RightDockWidgetArea, layer_control_panel_);
	ui_.menuView->addAction(layer_control_panel_->toggleViewAction());
	connect(rendering_layer_model_, SIGNAL(LayerPropertyChanged()), main_view_, SLOT(update()));

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

	/*vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName("./TestData/simple001.vtk");
	reader->SetReadAllScalars(1);
	reader->SetReadAllVectors(1);
	reader->Update();
	vtkDataObject* data = reader->GetOutput();
	vtkUnstructuredGrid* scatter_point_data = vtkUnstructuredGrid::SafeDownCast(data);
	vtkPointData* pData = scatter_point_data->GetPointData();
	vtkFloatArray* varray = vtkFloatArray::SafeDownCast(pData->GetArray(0));
	vtkPoints* p = scatter_point_data->GetPoints();

	int point_num = p->GetNumberOfPoints();
	int numArray = pData->GetNumberOfArrays();
	dataset_->original_point_values.resize(point_num);
	for (int i = 0; i < point_num; ++i) dataset_->original_point_values[i].resize(1);
	dataset_->original_point_pos.resize(point_num);
	for (int i = 0; i < point_num; ++i) dataset_->original_point_pos[i].resize(2);

	for (int i = 0; i < point_num; ++i){
		dataset_->original_point_pos[i][0] = p->GetPoint(i)[0];
		dataset_->original_point_pos[i][1] = p->GetPoint(i)[1];
		dataset_->original_point_values[i][0] = varray->GetValue(i);
	}*/

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
	/*vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName("./TestData/timestep0_0003600.vtk");
	reader->SetReadAllScalars(1);
	reader->SetReadAllVectors(1);
	reader->Update();
	vtkDataObject* data = reader->GetOutput();
	scatter_point_data_ = vtkUnstructuredGrid::SafeDownCast(data);

	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer< vtkPolyData >::New();
	polydata->SetPoints(scatter_point_data_->GetPoints());
	vtkPoints* points = vtkPoints::New();
	polydata->SetPoints(points);
	vtkCellArray* vert_array = vtkCellArray::New();
	polydata->SetVerts(vert_array);
	vtkIntArray* x_var_array = vtkIntArray::New();
	x_var_array->SetNumberOfComponents(1);
	vtkFloatArray* y_var_array = vtkFloatArray::New();
	y_var_array->SetNumberOfComponents(1);
	vtkFloatArray* xp_var_array = vtkFloatArray::New();
	xp_var_array->SetNumberOfComponents(1);
	vtkFloatArray* yp_var_array = vtkFloatArray::New();
	yp_var_array->SetNumberOfComponents(1);

	vtkIntArray* xarray = vtkIntArray::SafeDownCast(scatter_point_data_->GetPointData()->GetArray(0));
	vtkFloatArray* yarray = vtkFloatArray::SafeDownCast(scatter_point_data_->GetPointData()->GetArray(1));
	vtkFloatArray* xparray = vtkFloatArray::SafeDownCast(scatter_point_data_->GetPointData()->GetArray(2));
	vtkFloatArray* yparray = vtkFloatArray::SafeDownCast(scatter_point_data_->GetPointData()->GetArray(3));

	float x_scale1 = 0.45, x_scale2 = 0.48;
	float y_scale1 = 0.75, y_scale2 = 0.77;
	for (int i = 0; i < scatter_point_data_->GetPoints()->GetNumberOfPoints(); ++i) {
		float x = scatter_point_data_->GetPoint(i)[0];
		float y = scatter_point_data_->GetPoint(i)[1];
		double bounds[6];
		scatter_point_data_->GetBounds(bounds);
		if (x > bounds[0] + (bounds[1] - bounds[0]) * x_scale1 && x < bounds[0] + (bounds[1] - bounds[0]) * x_scale2
			&& y > bounds[2] + (bounds[3] - bounds[2]) * y_scale1 && y < bounds[2] + (bounds[3] - bounds[2]) * y_scale2) {
			points->InsertNextPoint(x, y, 0);
			vtkIdType id = points->GetNumberOfPoints() - 1;
			polydata->InsertNextCell(VTK_VERTEX, 1, &id);
			x_var_array->InsertNextTuple1(xarray->GetValue(i));
			y_var_array->InsertNextTuple1(yarray->GetValue(i));
			xp_var_array->InsertNextTuple1(xparray->GetValue(i));
			yp_var_array->InsertNextTuple1(yparray->GetValue(i));
		}
	}

	vtkSmartPointer< vtkUnstructuredGridWriter > writer = vtkSmartPointer< vtkUnstructuredGridWriter >::New();
	vtkUnstructuredGrid* ungrid = vtkUnstructuredGrid::New();
	ungrid->SetPoints(polydata->GetPoints());
	ungrid->GetPointData()->AddArray(y_var_array);
	int narray = ungrid->GetPointData()->GetNumberOfArrays();
	writer->SetInputData(ungrid);
	writer->SetFileName("./TestData/simple001.vtk");
	writer->Update();*/
}

void ScatterPointGlyph::OnActionCloseTriggered() {

}

void ScatterPointGlyph::OnActionExitTriggered() {
	
}

void ScatterPointGlyph::OnGlyphSelected(int x, int y) {
	double point[4];
	vtkInteractorObserver::ComputeDisplayToWorld(this->main_renderer_, (double)x, (double)this->main_view_->height() - y, 0, point);

	if (cluster_glyph_layer_ != NULL) {
		cluster_glyph_layer_->Select(point[0], point[1]);
		std::vector< bool > is_cluster_selected = cluster_glyph_layer_->GetClusterSelection();
		std::vector< int > cluster_color = cluster_point_rendering_layer_->GetClusterColor();
		parallel_dataset_->is_record_selected.resize(1);
		parallel_dataset_->is_record_selected[0].resize(dataset_->original_point_pos.size());
		parallel_dataset_->is_record_selected[0].assign(dataset_->original_point_pos.size(), false);
		parallel_dataset_->record_color.resize(1);
		parallel_dataset_->record_color[0].resize(dataset_->original_point_pos.size());

		for (int i = 0; i < is_cluster_selected.size(); ++i)
			if (is_cluster_selected[i]) {
				for (int j = 0; j < cluster_index.size(); ++j)
					if (cluster_index[j] == i) {
						parallel_dataset_->is_record_selected[0][j] = true;
						parallel_dataset_->record_color[0][j] = QColor(cluster_color[3 * i], cluster_color[3 * i + 1], cluster_color[3 * i + 2]);
					}
			}
	}
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

		HierarchicalTree* hier_tree = dynamic_cast<HierarchicalTree*>(cluster_tree_vec_[PERCEPTION_MODE]);
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
	std::vector< float > represent_value;
	vtkSmartPointer< vtkPoints > original_points = vtkSmartPointer< vtkPoints >::New();
	vtkSmartPointer< vtkPolyData > original_point_poly = vtkSmartPointer< vtkPolyData >::New();
	original_point_poly->SetPoints(original_points);
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i){
		double new_pos[3];
		new_pos[0] = dataset_->original_point_pos[i][0];
		new_pos[1] = dataset_->original_point_pos[i][1];
		new_pos[2] = 0;
		original_points->InsertNextPoint(new_pos);
		represent_value.push_back(dataset_->point_values[i][0]);
	}

	vtkSmartPointer< vtkVertexGlyphFilter > original_vertex_filter = vtkSmartPointer< vtkVertexGlyphFilter >::New();
	original_vertex_filter->SetInputData(original_point_poly);
	original_vertex_filter->Update();
	vtkSmartPointer< vtkPolyData > original_vertex_data = vtkSmartPointer< vtkPolyData >::New();
	original_vertex_data->ShallowCopy(original_vertex_filter->GetOutput());
	if (original_point_rendering_layer_ == NULL) {
		original_point_rendering_layer_ = new PointRenderingLayer;
		main_renderer_->AddActor(original_point_rendering_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Original Point Layer"), original_point_rendering_layer_->actor());
	}
	original_point_rendering_layer_->SetData(original_vertex_data);
	original_point_rendering_layer_->SetPointValue(represent_value);

	vtkSmartPointer< vtkPolyData > cluster_vert_data = vtkSmartPointer< vtkPolyData >::New();
	cluster_vert_data->DeepCopy(original_vertex_data);
	if (cluster_point_rendering_layer_ == NULL) {
		cluster_point_rendering_layer_ = new PointRenderingLayer;
		main_renderer_->AddActor(cluster_point_rendering_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Cluster Point Layer"), cluster_point_rendering_layer_->actor(), false);
	}
	cluster_point_rendering_layer_->SetData(cluster_vert_data);
	cluster_point_rendering_layer_->SetPointValue(represent_value);

	vtkSmartPointer< vtkPolyData > un_vert_data = vtkSmartPointer< vtkPolyData >::New();
	un_vert_data->DeepCopy(original_vertex_data);
	if (un_rendering_layer_ == NULL) {
		un_rendering_layer_ = new PointRenderingLayer;
		main_renderer_->AddActor(un_rendering_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Uncertainty Map"), un_rendering_layer_->actor(), false);
	}
	un_rendering_layer_->SetData(un_vert_data);

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

	main_renderer_->ResetCamera();
	main_view_->update();
}

void ScatterPointGlyph::OnClusterFinished() {
	if (cluster_glyph_layer_ == NULL) {
		cluster_glyph_layer_ = new ClusterGlyphLayer;
		main_renderer_->AddActor(cluster_glyph_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Hierarchical Cluster Glyph"), cluster_glyph_layer_->actor());
		cluster_glyph_layer_->SetData(dataset_);
	}

	if (sys_mode_ == PERCEPTION_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		TreeCommon* tree = cluster_tree_vec_[PERCEPTION_MODE];
		tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluter_num, cluster_index);
		cluster_glyph_layer_->SetRadiusRange(dis_per_pixel_ * 20, dis_per_pixel_ * 20);
		cluster_glyph_layer_->SetClusterIndex(cluter_num, cluster_index);

		cluster_point_rendering_layer_->SetClusterIndex(cluter_num, cluster_index);
	}

	if (sys_mode_ == IMMEDIATE_PERCEPTION_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		TreeCommon* tree = cluster_tree_vec_[IMMEDIATE_PERCEPTION_MODE];
		tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluter_num, cluster_index);
		cluster_glyph_layer_->SetRadiusRange(dis_per_pixel_ * 20, dis_per_pixel_ * 20);
		cluster_glyph_layer_->SetClusterIndex(cluter_num, cluster_index);

		cluster_point_rendering_layer_->SetClusterIndex(cluter_num, cluster_index);
	}

	if (sys_mode_ == IMMEDIATE_GESTALT_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		TreeCommon* tree = cluster_tree_vec_[IMMEDIATE_GESTALT_MODE];
		tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluter_num, cluster_index);
		cluster_glyph_layer_->SetRadiusRange(dis_per_pixel_ * 20, dis_per_pixel_ * 20);
		cluster_glyph_layer_->SetClusterIndex(cluter_num, cluster_index);

		cluster_point_rendering_layer_->SetClusterIndex(cluter_num, cluster_index);
	}

	if (sys_mode_ == UNCERTAINTY_MODE && cluster_tree_vec_[sys_mode_] != NULL) {
		UncertaintyTree* un_tree = dynamic_cast< UncertaintyTree* >(cluster_tree_vec_[UNCERTAINTY_MODE]);
		un_tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluter_num, cluster_index);
		cluster_glyph_layer_->SetRadiusRange(dis_per_pixel_ * 30, dis_per_pixel_ * 30);
		cluster_glyph_layer_->SetClusterIndex(cluter_num, cluster_index);

		cluster_point_rendering_layer_->SetClusterIndex(cluter_num, cluster_index);
		un_rendering_layer_->SetPointValue(un_tree->GetUncertainty());

		TransMapData* data = new TransMapData;
		data->node_center.resize(cluter_num);
		data->node_average_value.resize(cluter_num);
		data->node_num = cluter_num;
		data->var_num = dataset_->weights.size();
		for (int i = 0; i < cluter_num; ++i) {
			data->node_average_value[i].resize(data->var_num, 0);
			data->node_center[i].resize(2, 0);
		}

		std::vector< int > cluster_point_count;
		cluster_point_count.resize(cluter_num);
		cluster_point_count.assign(cluter_num, 0);
		for (int i = 0; i < cluster_index.size(); ++i) {
			int temp_index = cluster_index[i];
			if (temp_index < 0) continue;
			data->node_center[temp_index][0] += dataset_->original_point_pos[i][0];
			data->node_center[temp_index][1] += dataset_->original_point_pos[i][1];
			for (int j = 0; j < dataset_->weights.size(); ++j)
				data->node_average_value[temp_index][j] += (dataset_->original_point_values[i][j] - dataset_->original_value_ranges[j][0]) / (dataset_->original_value_ranges[j][1] - dataset_->original_value_ranges[j][0]);
			cluster_point_count[temp_index]++;
		}
		for (int i = 0; i < cluter_num; ++i)
			if (cluster_point_count[i] != 0) {
				data->node_center[i][0] /= cluster_point_count[i];
				data->node_center[i][1] /= cluster_point_count[i];
				for (int j = 0; j < dataset_->weights.size(); ++j)
					data->node_average_value[i][j] /= cluster_point_count[i];
			}
			else {
				std::cout << "Empty cluster detected." << std::endl;
				return;
			}

			vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
			for (int i = 0; i < data->node_num; ++i) {
				points->InsertNextPoint(data->node_center[i][0], data->node_center[i][1], 0);
			}
			vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
			polydata->SetPoints(points);
			vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer< vtkDelaunay2D >::New();
			delaunay->SetInputData(polydata);
			delaunay->SetTolerance(0.00001);
			delaunay->SetBoundingTriangulation(false);
			delaunay->Update();
			vtkPolyData* triangle_out = delaunay->GetOutput();

			data->node_connecting_status.resize(data->node_num);
			for (int k = 0; k < data->node_num; ++k) data->node_connecting_status[k].resize(data->node_num, false);
			for (int k = 0; k < triangle_out->GetNumberOfPolys(); ++k){
				vtkCell* cell = triangle_out->GetCell(k);
				int id1 = cell->GetPointId(0);
				int id2 = cell->GetPointId(1);
				int id3 = cell->GetPointId(2);
				data->node_connecting_status[id1][id2] = true;
				data->node_connecting_status[id2][id1] = true;
				data->node_connecting_status[id2][id3] = true;
				data->node_connecting_status[id3][id2] = true;
				data->node_connecting_status[id1][id3] = true;
				data->node_connecting_status[id3][id1] = true;
			}

			data->var_repsentative_color = cluster_point_rendering_layer_->GetClusterColor();

			trans_map_->SetNodeRadius(dis_per_pixel_ * 20);
			trans_map_->SetData(data);
			trans_map_->SetEnabled(true);
	}

	main_view_->update();
}

float ScatterPointGlyph::GetMainViewDisPerPixel() {
	int window_width = this->main_view_->width();
	vtkCamera* camera = this->main_renderer_->GetActiveCamera();
	vtkMatrix4x4* mMatrix = camera->GetViewTransformMatrix();
	vtkMatrix4x4* pMatrix = camera->GetProjectionTransformMatrix(this->main_renderer_);
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