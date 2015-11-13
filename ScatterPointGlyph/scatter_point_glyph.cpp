#include "scatter_point_glyph.h"

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
#include "gestalt_processor.h"

//#define TEST_WRITE

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), scatter_point_data_(NULL), sys_mode_(PERCEPTION_MODE), hier_solver_(NULL),
	cluster_glyph_layer_(NULL), point_rendering_layer_(NULL), original_rendering_layer_(NULL),
		gestalt_processor_(NULL) {
	ui_.setupUi(this);

	this->InitWidget();
}

ScatterPointGlyph::~ScatterPointGlyph() {

}

void ScatterPointGlyph::InitWidget() {
	main_view_ = new QVTKWidget;
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
	connect(ui_.action_perception_driven, SIGNAL(triggered()), this, SLOT(OnActionPerceptionDrivenTriggered()));
	connect(hier_para_widget_, SIGNAL(ClusterNumberChanged(int)), this, SLOT(OnHierClusterNumberChanged(int)));
}

void ScatterPointGlyph::OnActionOpenTriggered() {
	if (scatter_point_data_ != NULL) {
		QMessageBox::warning(this, tr("Warning"), tr("Data has been opened!"));
		return;
	}

	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
#ifdef TEST_WRITE
	reader->SetFileName("./TestData/timestep0_0003600.vtk");
#else
	reader->SetFileName("./TestData/simple001.vtk");
#endif
	reader->SetReadAllScalars(1);
	reader->SetReadAllVectors(1);
	reader->Update();
	vtkDataObject* data = reader->GetOutput();
	scatter_point_data_ = vtkUnstructuredGrid::SafeDownCast(data);
	vtkPointData* pData = scatter_point_data_->GetPointData();
	int numArray = pData->GetNumberOfArrays();

#ifdef TEST_WRITE
	vtkIntArray* xarray = vtkIntArray::SafeDownCast(pData->GetArray(0));
	std::cout << pData->GetArray(0)->GetName() << std::endl;
	vtkFloatArray* varray = vtkFloatArray::SafeDownCast(pData->GetArray(1));
	std::cout << pData->GetArray(1)->GetName() << std::endl;
	vtkFloatArray* xparray = vtkFloatArray::SafeDownCast(pData->GetArray(2));
	std::cout << pData->GetArray(2)->GetName() << std::endl;
	vtkFloatArray* yparray = vtkFloatArray::SafeDownCast(pData->GetArray(3));
	std::cout << pData->GetArray(3)->GetName() << std::endl;
#else
	vtkFloatArray* varray = vtkFloatArray::SafeDownCast(pData->GetArray(0));
#endif
	vtkPoints* p = scatter_point_data_->GetPoints();

	layer_control_widget_->SetDatasetInfo(QString("Point"), p->GetNumberOfPoints());

	/// construct the data
	int sample_num = p->GetNumberOfPoints();
	variable_weight_.resize(numArray, 1.0 / numArray);
	point_values_.resize(sample_num);
	for (int i = 0; i < sample_num; ++i) point_values_[i].resize(1);
	point_pos_.resize(sample_num);
	for (int i = 0; i < sample_num; ++i) point_pos_[i].resize(2);

	for (int i = 0; i < sample_num; ++i){
		point_pos_[i][0] = p->GetPoint(i)[0];
		point_pos_[i][1] = p->GetPoint(i)[1];
		point_values_[i][0] = varray->GetValue(i);
	}

	NormalizePosition(point_pos_);
	NormalizeVector(point_values_);

	this->AddPointData2View();

#ifndef TEST_WRITE
	this->PreProcess();
#endif
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
		connect(hier_solver_, SIGNAL(finished()), this, SLOT(OnOneStepHierFinished()));
		connect(hier_solver_, SIGNAL(CombinedClusterChanged(int, int)), this, SLOT(OnCombinedClusterUpdated(int, int)));
	}

	/*ContinuityExtractor* con_extractor = new ContinuityExtractor;
	con_extractor->SetData(x, y);

	hier_solver_->SetData(data_value, weight);
	expected_cluster_num_ = 200;
	hier_para_widget_->SetMaxClusterNumber(expected_cluster_num_);
	hier_solver_->SetInitialClusterNum(expected_cluster_num_);
	this->OnOneStepHierFinished();*/
}

void ScatterPointGlyph::PerceptionPreProcess() {
	if (gestalt_processor_ == NULL) {
		gestalt_processor_ = new GestaltProcessor;
	}

	gestalt_processor_->SetData(point_pos_, point_values_, variable_weight_);

#ifdef DEBUG_ON
	float min_dis = 1e10;
	for (int i = 0; i < point_pos_.size() - 1; ++i)
		for (int j = i + 1; j < point_pos_.size(); ++j) {
			float temp_dis = abs(point_pos_[i][0] - point_pos_[j][0]) + abs(point_pos_[i][1] - point_pos_[j][1]);
			if (temp_dis < min_dis) min_dis = temp_dis;
		}
#endif

	connect(gestalt_processor_, SIGNAL(FinalGlyphUpdated()), this, SLOT(OnGestaltUpdated()));
	connect(gestalt_processor_, SIGNAL(KMeansClusterFinished()), this, SLOT(OnKmeansClusterFinished()));

	/*std::vector< float > threshold;
	threshold.push_back(0.4);
	threshold.push_back(0.01);
	gestalt_processor_->SetGestaltThreshold(threshold);
	gestalt_processor_->SetDistanceThreshold(0.02);
	gestalt_processor_->GenerateLayout();*/
}

void ScatterPointGlyph::ExecPerceptionClustering() {
	
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

}

void ScatterPointGlyph::OnActionPerceptionDrivenTriggered() {
	int window_width = this->main_view_->width();
	vtkCamera* camera = this->main_renderer_->GetActiveCamera();
	
	vtkMatrix4x4* mMatrix = camera->GetViewTransformMatrix();
	vtkMatrix4x4* pMatrix = camera->GetProjectionTransformMatrix(this->main_renderer_);
	double point_one[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 1, 0, 0, point_one);
	double point_two[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 2, 0, 0, point_two);


	float width_per_scale = abs(point_one[0] - point_two[0]);
	std::vector< float > threshold;
	threshold.push_back(0.4);
	threshold.push_back(30 / width_per_scale);
	gestalt_processor_->SetGestaltThreshold(threshold);
	gestalt_processor_->SetDistanceThreshold(60 / width_per_scale);
	gestalt_processor_->GenerateLayout();

	main_view_->update();
} 

void ScatterPointGlyph::AddPointData2View() {
	vtkSmartPointer< vtkActor > test_actor = vtkSmartPointer< vtkActor >::New();
	vtkSmartPointer< vtkPolyDataMapper > test_mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
	vtkSmartPointer< vtkPolyData > test_poly = vtkSmartPointer< vtkPolyData >::New();

#ifdef TEST_WRITE
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
#else
	vtkSmartPointer< vtkPoints > new_points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < point_pos_.size(); ++i)
		new_points->InsertNextPoint(point_pos_[i][0], point_pos_[i][1], 0);
	test_poly->SetPoints(new_points);
	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
	vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(test_poly);
	vertexFilter->Update();
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->ShallowCopy(vertexFilter->GetOutput());
#endif

	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(1);
	colors->SetName("Colors");
	for (size_t i = 0; i < polydata->GetPoints()->GetNumberOfPoints(); ++i){
		colors->InsertNextTuple1(128);
	}
	polydata->GetPointData()->SetScalars(colors);

#ifdef TEST_WRITE
	vtkSmartPointer< vtkUnstructuredGridWriter > writer = vtkSmartPointer< vtkUnstructuredGridWriter >::New();
	vtkUnstructuredGrid* ungrid = vtkUnstructuredGrid::New();
	ungrid->SetPoints(polydata->GetPoints());
	ungrid->GetPointData()->AddArray(y_var_array);
	int narray = ungrid->GetPointData()->GetNumberOfArrays();
	writer->SetInputData(ungrid);
	writer->SetFileName("./TestData/simple001.vtk");
	writer->Update();
#endif

	if (point_rendering_layer_ == NULL) point_rendering_layer_ = new PointRenderingLayer;
	point_rendering_layer_->SetData(polydata);

	if (original_rendering_layer_ == NULL) original_rendering_layer_ = new PointRenderingLayer;
	vtkPolyData* data_copy = vtkPolyData::New();
	data_copy->DeepCopy(polydata);
	for (int i = 0; i < data_copy->GetPoints()->GetNumberOfPoints(); ++i){
		double* pos = data_copy->GetPoint(i);
		double new_pos[3];
		new_pos[0] = pos[0];
		new_pos[1] = pos[1] + 1;
		new_pos[2] = pos[2];
		data_copy->GetPoints()->SetPoint(i, new_pos);
	}
	original_rendering_layer_->SetData(data_copy);
	main_renderer_->AddActor(original_rendering_layer_->actor());
	rendering_layer_model_->AddLayer(QString("Original Point Layer"), original_rendering_layer_->actor());
	std::vector< float > temp_value;
	for (int i = 0; i < point_values_.size(); ++i) temp_value.push_back(point_values_[i][0]);
	original_rendering_layer_->SetPointValue(temp_value);

	rendering_layer_model_->AddLayer(QString("Point Layer"), point_rendering_layer_->actor());

	main_renderer_->AddActor(point_rendering_layer_->actor());
	main_renderer_->ResetCamera();

	if (cluster_glyph_layer_ == NULL) {
		cluster_glyph_layer_ = new ClusterGlyphLayer;
		main_renderer_->AddActor(cluster_glyph_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Hier Cluster Glyph"), cluster_glyph_layer_->actor());
	}

	cluster_glyph_layer_->SetData(point_pos_, point_values_);

	main_view_->update();
}

void ScatterPointGlyph::OnHierClusterNumberChanged(int num) {
	expected_cluster_num_ = num;
	hier_solver_->SetClusterNum(num);

	if (hier_solver_->GetCurrentClusterNum() > expected_cluster_num_)
		hier_solver_->start();
}

void ScatterPointGlyph::OnCombinedClusterUpdated(int cluster_one, int cluster_two) {
	cluster_glyph_layer_->SetHighlighCluster(cluster_one, cluster_two);
}

void ScatterPointGlyph::OnOneStepHierFinished() {
	/// TODO: update the visualization
}

void ScatterPointGlyph::OnGestaltUpdated() {
	std::vector< int > point_index;
	gestalt_processor_->GetFinalGlyphPoint(point_index);
	// update the point layer
	point_rendering_layer_->SetHighlightPointIndex(point_index);

	// add glyph to the glyph layer
	//cluster_glyph_layer_->AddGestaltGlyph(point_index);
}

void ScatterPointGlyph::OnKmeansClusterFinished() {
	std::vector< std::vector< int > > index;
	gestalt_processor_->GetKmeansClusterResult(index);
	for (int i = 0; i < index.size(); ++i) original_rendering_layer_->SetHighlightPointIndex(index[i]);
}

void ScatterPointGlyph::NormalizeVector(std::vector< std::vector< float > >& vec){
	
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
			vec[j][i] = (vec[j][i] - ranges[i][0]) / max_range;
	}

}