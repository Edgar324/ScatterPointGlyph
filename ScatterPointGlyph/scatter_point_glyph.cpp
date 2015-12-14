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

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), dataset_(NULL), sys_mode_(PERCEPTION_MODE), 
	cluster_glyph_layer_(NULL), original_point_rendering_layer_(NULL), cluster_point_rendering_layer_(NULL),
	map_rendering_layer_(NULL), dis_per_pixel_(0.0), min_pixel_radius_(5) {

	ui_.setupUi(this);

	data_manager_ = WrfDataManager::GetInstance();

	this->InitWidget();

	cluster_tree_vec_.resize(2, NULL);
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
	main_renderer_->SetBackground(1.0, 1.0, 1.0);
	main_view_->GetRenderWindow()->AddRenderer(main_renderer_);
	connect(main_view_, SIGNAL(ViewUpdated()), this, SLOT(UpdateClusterView()));

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
	sys_mode_action_group_->setExclusive(true);

	connect(ui_.actionVTK_Unstructured_Grid_Data, SIGNAL(triggered()), this, SLOT(OnActionOpenVtkFileTriggered()));
	connect(ui_.actionRaw_Grid_Data, SIGNAL(triggered()), this, SLOT(OnActionOpenRawGridFileTriggered()));
	connect(ui_.actionExec, SIGNAL(triggered()), this, SLOT(UpdateClusterView()));
}

void ScatterPointGlyph::OnActionOpenVtkFileTriggered() {
	if (dataset_ == NULL) dataset_ = new ScatterPointDataset;

	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
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
	}

	dataset_->DirectConstruct();
	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionOpenRawGridFileTriggered() {
	if (dataset_ == NULL) dataset_ = new ScatterPointDataset;

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

void ScatterPointGlyph::OnSysmodeChanged() {
	if (ui_.action_perception_driven->isChecked()) {
		sys_mode_ = PERCEPTION_MODE;
	} else if (ui_.action_hierarchical_clustering->isChecked()) {
		sys_mode_ = HIER_MODE;
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
	default:
		break;
	}
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
		rendering_layer_model_->AddLayer(QString("Original Point Layer"), original_point_rendering_layer_->actor(), false);
	}
	original_point_rendering_layer_->SetData(original_vertex_data);
	original_point_rendering_layer_->SetPointValue(dataset_->point_values);

	vtkSmartPointer< vtkPolyData > cluster_vert_data = vtkSmartPointer< vtkPolyData >::New();
	cluster_vert_data->DeepCopy(original_vertex_data);
	if (cluster_point_rendering_layer_ == NULL) {
		cluster_point_rendering_layer_ = new PointRenderingLayer;
		main_renderer_->AddActor(cluster_point_rendering_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Cluster Point Layer"), cluster_point_rendering_layer_->actor());
	}
	cluster_point_rendering_layer_->SetData(cluster_vert_data);
	cluster_point_rendering_layer_->SetPointValue(dataset_->point_values);

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
		int cluter_num;
		std::vector< int > cluster_index;

		TreeCommon* tree = cluster_tree_vec_[PERCEPTION_MODE];
		tree->GetClusterResult(this->dis_per_pixel_ / (dataset_->original_pos_ranges[0][1] - dataset_->original_pos_ranges[0][0]), cluter_num, cluster_index);
		cluster_glyph_layer_->SetRadiusRange(dis_per_pixel_ * 20, dis_per_pixel_ * 20);
		cluster_glyph_layer_->SetClusterIndex(cluter_num, cluster_index);

		cluster_point_rendering_layer_->SetClusterIndex(cluter_num, cluster_index);
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