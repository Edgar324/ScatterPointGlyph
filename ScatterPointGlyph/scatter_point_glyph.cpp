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

#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDockWidget>

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

//#define SAVE_TXT_FILE
//#define WEATHREVIS

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), dataset_(NULL), sys_mode_(PERCEPTION_MODE), expected_cluster_num_(30),
	cluster_glyph_layer_(NULL), original_point_rendering_layer_(NULL), sample_point_rendering_layer_(NULL),
	map_rendering_layer_(NULL),	gestalt_processor_(NULL), hier_solver_(NULL), current_solver_(NULL) {

	ui_.setupUi(this);

	data_manager_ = WrfDataManager::GetInstance();

	this->InitWidget();
}

ScatterPointGlyph::~ScatterPointGlyph() {

}

void ScatterPointGlyph::InitWidget() {
	main_view_ = new ScatterPointView;
	main_view_->setFocusPolicy(Qt::StrongFocus);
	QHBoxLayout* central_widget_layout = new QHBoxLayout;
	central_widget_layout->addWidget(main_view_);
	ui_.centralWidget->setLayout(central_widget_layout);

	main_renderer_ = vtkRenderer::New();
	main_renderer_->RemoveAllViewProps();
	main_view_->GetRenderWindow()->AddRenderer(main_renderer_);
	main_renderer_->SetBackground(1.0, 1.0, 1.0);

	layer_control_widget_ = new LayerControlWidget;
	rendering_layer_model_ = new RenderingLayerModel;
	layer_control_widget_->SetRenderingLayerModel(rendering_layer_model_);
	layer_control_panel_ = new QDockWidget(QString("Layer Control Panel"), this);
	layer_control_panel_->setWidget(layer_control_widget_);
	this->addDockWidget(Qt::RightDockWidgetArea, layer_control_panel_);
	ui_.menuView->addAction(layer_control_panel_->toggleViewAction());

	hier_para_widget_ = new HierParaWidget;
	hier_para_panel_ = new QDockWidget(QString("Hierarchical Clustering Panel"), this);
	hier_para_panel_->setWidget(hier_para_widget_);
	this->addDockWidget(Qt::RightDockWidgetArea, hier_para_panel_);
	ui_.menuView->addAction(hier_para_panel_->toggleViewAction());

	connect(ui_.action_open, SIGNAL(triggered()), this, SLOT(OnActionOpenTriggered()));
	connect(ui_.action_hierarchical_clustering, SIGNAL(triggered()), this, SLOT(OnActionHierarchicalClusteringTriggered()));
	connect(ui_.action_perception_driven, SIGNAL(triggered()), this, SLOT(OnActionPerceptionDrivenTriggered()));
	connect(hier_para_widget_, SIGNAL(ClusterNumberChanged(int)), this, SLOT(OnHierClusterNumberChanged(int)));
	connect(main_view_, SIGNAL(ViewUpdated()), this, SLOT(OnMainViewUpdated()));
}

void ScatterPointGlyph::OnActionOpenTriggered() {
	if (scatter_point_data_ != NULL) {
		QMessageBox::warning(this, tr("Warning"), tr("Data has been opened!"));
		return;
	}

	if (dataset_ == NULL) dataset_ = new ScatterPointDataset;
	dataset_->weights.resize(1, 1.0);

#ifdef WEATHREVIS
	data_manager_->LoadEnsembleData(WRF_ACCUMULATED_PRECIPITATION, std::string("E:/Data/ens_l/apcp_sfc_latlon_all_20150401_20150430_liaoLeeVC4.nc"));
	std::vector< WrfGridValueMap* > maps;
	QDateTime temp_datetime = QDateTime(QDate(2015, 4, 5), QTime(0, 0));
	data_manager_->GetGridValueMap(temp_datetime, WRF_NCEP_ENSEMBLES, WRF_ACCUMULATED_PRECIPITATION, 24, maps);
	MapRange map_range;
	data_manager_->GetModelMapRange(WRF_NCEP_ENSEMBLES, map_range);

	int point_num = map_range.x_grid_number * map_range.y_grid_number;
	dataset_->original_point_values.resize(point_num);
	for (int i = 0; i < point_num; ++i) dataset_->original_point_values[i].resize(1);
	dataset_->original_point_pos.resize(point_num);
	for (int i = 0; i < point_num; ++i) dataset_->original_point_pos[i].resize(2);

	for (int i = 0; i < point_num; ++i){
		int w = i % map_range.x_grid_number;
		int h = i / map_range.x_grid_number;
		dataset_->original_point_pos[i][0] = map_range.start_x + map_range.x_grid_space * w;
		dataset_->original_point_pos[i][1] = map_range.start_y + map_range.y_grid_space * h;
		dataset_->original_point_values[i][0] = maps[0]->values[i];
	}
	dataset_->is_structured_data = true;
	dataset_->w = map_range.x_grid_number;
	dataset_->h = map_range.y_grid_number;

	NormalizeValues(dataset_->original_point_values);
	dataset_->DirectConstruct();
	NormalizePosition(dataset_->point_pos);

	/*if (map_rendering_layer_ == NULL) {
		map_rendering_layer_ = new MapRenderingLayer;
		map_rendering_layer_->SetFieldData(maps[0]->values, map_range.x_grid_number, map_range.y_grid_number, map_range.x_grid_space, map_range.y_grid_space, map_range.start_x, map_range.start_y);
		main_renderer_->AddActor(map_rendering_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Map"), map_rendering_layer_->actor());
	}*/
#else
	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName("./TestData/simple001.vtk");
	reader->SetReadAllScalars(1);
	reader->SetReadAllVectors(1);
	reader->Update();
	vtkDataObject* data = reader->GetOutput();
	scatter_point_data_ = vtkUnstructuredGrid::SafeDownCast(data);
	vtkPointData* pData = scatter_point_data_->GetPointData();
	vtkFloatArray* varray = vtkFloatArray::SafeDownCast(pData->GetArray(0));
	vtkPoints* p = scatter_point_data_->GetPoints();

	layer_control_widget_->SetDatasetInfo(QString("Point"), p->GetNumberOfPoints());

	int point_num = p->GetNumberOfPoints();
	int numArray = pData->GetNumberOfArrays();
	dataset_->original_point_values.resize(point_num - 1);
	for (int i = 0; i < point_num - 1; ++i) dataset_->original_point_values[i].resize(1);
	dataset_->original_point_pos.resize(point_num - 1);
	for (int i = 0; i < point_num - 1; ++i) dataset_->original_point_pos[i].resize(2);

	for (int i = 0; i < point_num; ++i){
		if (i == 31) continue;
		if (i > 31) {
			dataset_->original_point_pos[i - 1][0] = p->GetPoint(i)[0];
			dataset_->original_point_pos[i - 1][1] = p->GetPoint(i)[1];
			dataset_->original_point_values[i - 1][0] = varray->GetValue(i);
		}
		else {
			dataset_->original_point_pos[i][0] = p->GetPoint(i)[0];
			dataset_->original_point_pos[i][1] = p->GetPoint(i)[1];
			dataset_->original_point_values[i][0] = varray->GetValue(i);
		}
	}

	NormalizePosition(dataset_->original_point_pos);
	NormalizeValues(dataset_->original_point_values);
	dataset_->DirectConstruct();
#endif

#ifdef SAVE_TXT_FILE
	std::ofstream output("E:/GeoVis/scatterdatareview/1.txt");
	if (output.good()) {
		output << point_num - 1 << std::endl;
		for (int i = 0; i < point_num - 1; ++i) {
			output << dataset_->original_point_pos[i][0] << std::endl;
			output << dataset_->original_point_pos[i][1] << std::endl;
			output << dataset_->original_point_values[i][0] << std::endl;
		}
	}
	output.close();
#endif

	this->AddPointData2View();
	this->PreProcess();
}

void ScatterPointGlyph::OnActionExtractDataTriggered() {
	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
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
	writer->Update();
}

void ScatterPointGlyph::PreProcess() {
	switch (sys_mode_)
	{
	case ScatterPointGlyph::HIER_MODE:
		HierarchicalPreProcess();
		break;
	case ScatterPointGlyph::PERCEPTION_MODE:
		PerceptionPreProcess();
		break;
	default:
		break;
	}
}

void ScatterPointGlyph::HierarchicalPreProcess() {
	if (hier_solver_ == NULL) {
		hier_solver_ = new HierSolver;
		connect(hier_solver_, SIGNAL(ClusterUpdated(int)), this, SLOT(OnClusterAggregated(int)));
		connect(hier_solver_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
	}

	hier_solver_->SetData(dataset_);
	hier_solver_->SetExpectedClusterNum(expected_cluster_num_);
}

void ScatterPointGlyph::PerceptionPreProcess() {
	if (gestalt_processor_ == NULL) {
		gestalt_processor_ = new GestaltProcessor2;
		connect(gestalt_processor_, SIGNAL(ClusterUpdated(int)), this, SLOT(OnClusterAggregated(int)));
		connect(gestalt_processor_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
	}
}

void ScatterPointGlyph::OnActionCloseTriggered() {
	if (scatter_point_data_ != NULL) {
		scatter_point_data_->Delete();
		scatter_point_data_ = NULL;
	}
	/// TODO: clear the memory
}

void ScatterPointGlyph::OnActionExitTriggered() {
	
}

void ScatterPointGlyph::OnActionHierarchicalClusteringTriggered() {
	current_solver_ = hier_solver_;
	hier_solver_->start();
}

void ScatterPointGlyph::OnActionPerceptionDrivenTriggered() {
	current_solver_ = gestalt_processor_;

	int window_width = this->main_view_->width();
	vtkCamera* camera = this->main_renderer_->GetActiveCamera();
	
	vtkMatrix4x4* mMatrix = camera->GetViewTransformMatrix();
	vtkMatrix4x4* pMatrix = camera->GetProjectionTransformMatrix(this->main_renderer_);
	double point_one[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 1, 0, 0, point_one);
	double point_two[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 2, 0, 0, point_two);

	float width_per_scale = abs(point_one[0] - point_two[0]);
	/*std::vector< float > threshold;
	threshold.push_back(0.4);
	threshold.push_back(30 / width_per_scale);
	gestalt_processor_->SetDisThreshold(20 / width_per_scale);
	cluster_glyph_layer_->SetRadiusRange(20 / width_per_scale, 20 / width_per_scale * 0.1);
	gestalt_processor_->start();*/
	HierarchicalTree* tree = new HierarchicalTree(dataset_);
	std::vector< float > weights;
	weights.push_back(1.0);
	tree->SetVariableWeights(weights);
	tree->GenerateCluster(1.0 / width_per_scale, 2);

	main_view_->update();
} 

void ScatterPointGlyph::AddPointData2View() {
	vtkSmartPointer< vtkPoints > original_points = vtkSmartPointer< vtkPoints >::New();
	vtkSmartPointer< vtkPolyData > original_point_poly = vtkSmartPointer< vtkPolyData >::New();
	original_point_poly->SetPoints(original_points);
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i){
		double new_pos[3];
		new_pos[0] = dataset_->original_point_pos[i][0];
		new_pos[1] = dataset_->original_point_pos[i][1];
		new_pos[2] = 0;
		original_points->InsertNextPoint(new_pos);
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
	//original_point_rendering_layer_->SetPointValue(dataset_->point_values);

	vtkSmartPointer< vtkPolyData > sample_data = vtkSmartPointer< vtkPolyData >::New();
	vtkSmartPointer< vtkPoints > sample_points = vtkSmartPointer< vtkPoints >::New();
	sample_data->SetPoints(sample_points);
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i){
		double new_pos[3];
		new_pos[0] = dataset_->original_point_pos[i][0];
		new_pos[1] = dataset_->original_point_pos[i][1] + 1;
		new_pos[2] = 0;
		sample_points->InsertNextPoint(new_pos);
	}
	vtkSmartPointer< vtkVertexGlyphFilter > sample_vertex_filter = vtkSmartPointer< vtkVertexGlyphFilter >::New();
	sample_vertex_filter->SetInputData(sample_data);
	sample_vertex_filter->Update();
	vtkSmartPointer< vtkPolyData > sample_vert_data = vtkSmartPointer< vtkPolyData >::New();
	sample_vert_data->ShallowCopy(sample_vertex_filter->GetOutput());

	if (sample_point_rendering_layer_ == NULL) {
		sample_point_rendering_layer_ = new PointRenderingLayer;
		main_renderer_->AddActor(sample_point_rendering_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Sample Point Layer"), sample_point_rendering_layer_->actor());
	}
	sample_point_rendering_layer_->SetData(sample_vert_data);
	//sample_point_rendering_layer_->SetPointValue(dataset_->point_values);

	main_renderer_->ResetCamera();

	if (cluster_glyph_layer_ == NULL) {
		cluster_glyph_layer_ = new ClusterGlyphLayer;
		main_renderer_->AddActor(cluster_glyph_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Hierarchical Cluster Glyph"), cluster_glyph_layer_->actor());
	}

	cluster_glyph_layer_->SetData(dataset_);

	main_view_->update();
}

void ScatterPointGlyph::OnMainViewUpdated() {
	if (current_solver_ == NULL) return;
	switch (sys_mode_)
	{
	case HIER_MODE:
		OnActionHierarchicalClusteringTriggered();
		break;
	case PERCEPTION_MODE:
		OnActionPerceptionDrivenTriggered();
		break;
	default:
		break;
	}
}

void ScatterPointGlyph::OnHierClusterNumberChanged(int num) {
	expected_cluster_num_ = num;
	hier_solver_->SetExpectedClusterNum(num);
	hier_solver_->start();
}

void ScatterPointGlyph::OnClusterAggregated(int cluster_index) {
	std::vector< int > point_index;
	current_solver_->GetFinalClusterPoint(cluster_index, point_index);

	// update the point layer
	sample_point_rendering_layer_->SetHighlightPointIndex(point_index);
}

void ScatterPointGlyph::OnClusterFinished() {
	int cluster_count = 0;
	std::vector< int > cluster_index;
	current_solver_->GetCluster(cluster_count, cluster_index);

	original_point_rendering_layer_->SetClusterIndex(cluster_count, cluster_index);
	// add glyph to the glyph layer
	//cluster_glyph_layer_->SetClusterIndex(cluster_count, cluster_index);

	main_view_->update();
}

void ScatterPointGlyph::NormalizeValues(std::vector< std::vector< float > >& vec){
	for (int i = 0; i < vec[0].size(); ++i){
		float minValue = 1e10;
		float maxValue = -1e10;

		for (int j = 0; j < vec.size(); ++j){
			if (minValue > vec[j][i]) minValue = vec[j][i];
			if (maxValue < vec[j][i]) maxValue = vec[j][i];
		}

		if (maxValue - minValue != 0) {
			for (int j = 0; j < vec.size(); ++j)
				vec[j][i] = (vec[j][i] - minValue) / (maxValue - minValue);
		}
		else {
			for (int j = 0; j < vec.size(); ++j) vec[j][i] = 0.5;
		}
	}
}

void ScatterPointGlyph::NormalizePosition(std::vector< std::vector< float > >& vec) {
	float max_range = -1e10;
	std::vector< std::vector< float > > ranges;
	ranges.resize(vec[0].size());
	for (int i = 0; i < vec[0].size(); ++i){
		float minValue = 1e10;
		float maxValue = -1e10;

		for (int j = 0; j < vec.size(); ++j){
			if (minValue > vec[j][i]) minValue = vec[j][i];
			if (maxValue < vec[j][i]) maxValue = vec[j][i];
		}
		if (maxValue - minValue > max_range) max_range = maxValue - minValue;

		ranges[i].push_back(minValue);
		ranges[i].push_back(maxValue);
	}
	for (int i = 0; i < vec[0].size(); ++i){
		for (int j = 0; j < vec.size(); ++j)
			vec[j][i] = (vec[j][i] - (ranges[i][0] + ranges[i][1]) / 2) / max_range + 0.5;
	}
}