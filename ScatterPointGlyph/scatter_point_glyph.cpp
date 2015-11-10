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

#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDockWidget>

#include "layer_control_widget.h"
#include "rendering_layer_model.h"
#include "hier_para_widget.h"
#include "hier_solver.h"
#include "hier_cluster_glyph_layer.h"
#include "continuity_extractor.h"

//#define TEST_WRITE

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), scatter_point_data_(NULL), sys_mode_(HIER_MODE), hier_solver_(NULL),
		hier_cluster_glyph_layer_(NULL) {
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

#ifndef TEST_WRITE
	this->PreProcess();
#endif
	this->AddPointData2View();
}

void ScatterPointGlyph::PreProcess() {
	switch (sys_mode_)
	{
	case ScatterPointGlyph::HIER_MODE:
		HierarchicalPreProcess();
		break;
	case ScatterPointGlyph::PERCEPTION_MODE:
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

	vtkPointData* pData = scatter_point_data_->GetPointData();
	int numArray = pData->GetNumberOfArrays();
	vtkIntArray* xarray = vtkIntArray::SafeDownCast(pData->GetArray(0));
	vtkFloatArray* yarray = vtkFloatArray::SafeDownCast(pData->GetArray(1));
	vtkFloatArray* xparray = vtkFloatArray::SafeDownCast(pData->GetArray(2));
	vtkFloatArray* yparray = vtkFloatArray::SafeDownCast(pData->GetArray(3));
	vtkPoints* p = scatter_point_data_->GetPoints();

	/// construct the data
	std::vector< std::vector< float > > data_value;
	std::vector< float > weight;

	int sample_num = p->GetNumberOfPoints();
	weight.resize(2, 1);
	data_value.resize(sample_num);
	for (int i = 0; i < sample_num; ++i) data_value[i].resize(2);

	std::vector< float > x;
	std::vector< float > y;
	x.resize(sample_num);
	y.resize(sample_num);

	for (int i = 0; i < sample_num; ++i){
		data_value[i][0] = p->GetPoint(i)[0];
		data_value[i][1] = p->GetPoint(i)[1];

		x[i] = p->GetPoint(i)[0];
		y[i] = p->GetPoint(i)[1];
		/*data_value[i][2] = xarray->GetValue(i);
		data_value[i][3] = yarray->GetValue(i);
		data_value[i][4] = xparray->GetValue(i);
		data_value[i][5] = yparray->GetValue(i);*/
	}

	ContinuityExtractor* con_extractor = new ContinuityExtractor;
	con_extractor->SetData(x, y);

	hier_solver_->SetData(data_value, weight);
	expected_cluster_num_ = 200;
	hier_para_widget_->SetMaxClusterNumber(expected_cluster_num_);
	hier_solver_->SetInitialClusterNum(expected_cluster_num_);
	this->OnOneStepHierFinished();

	layer_control_widget_->SetDatasetInfo(QString("Point"), p->GetNumberOfPoints());
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
		}
	}
#else
	test_poly->SetPoints(scatter_point_data_->GetPoints());
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
	vtkSmartPointer< vtkUnstructuredGrid > new_un_grid = vtkSmartPointer< vtkUnstructuredGrid >::New();
	new_un_grid->SetPoints(polydata->GetPoints());
	writer->SetInputData(new_un_grid);
	writer->SetFileName("./TestData/simple001.vtk");
	writer->Update();
#endif

	test_mapper->SetInputData(polydata);
	test_actor->SetMapper(test_mapper);
	test_actor->GetProperty()->SetPointSize(3);

	rendering_layer_model_->AddLayer(QString("Point Layer"), test_actor);

	main_renderer_->AddActor(test_actor);
	main_renderer_->ResetCamera();
	main_view_->update();
}

void ScatterPointGlyph::OnHierClusterNumberChanged(int num) {
	expected_cluster_num_ = num;
	hier_solver_->SetClusterNum(num);

	if (hier_solver_->GetCurrentClusterNum() > expected_cluster_num_)
		hier_solver_->start();
}

void ScatterPointGlyph::OnCombinedClusterUpdated(int cluster_one, int cluster_two) {
	hier_cluster_glyph_layer_->SetHighlighCluster(cluster_one, cluster_two);
	main_view_->update();
}

void ScatterPointGlyph::OnOneStepHierFinished() {
	/// TODO: update the visualization
	if (hier_cluster_glyph_layer_ == NULL) {
		hier_cluster_glyph_layer_ = new HierClusterGlyphLayer;
		main_renderer_->AddActor(hier_cluster_glyph_layer_->actor());
		rendering_layer_model_->AddLayer(QString("Hier Cluster Glyph"), hier_cluster_glyph_layer_->actor());
	}

	std::vector< float > pos;
	std::vector< std::vector< float > > value;
	hier_solver_->GetGlyphData(pos, value);
	hier_cluster_glyph_layer_->SetData(pos, value);
	main_view_->update();

	if (hier_solver_->GetCurrentClusterNum() > expected_cluster_num_) 
		hier_solver_->start();
}