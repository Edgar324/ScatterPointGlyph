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

#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QMessageBox>

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), scatter_point_data_(NULL) {
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
	context_menu_ = new QMenu;

	connect(ui_.action_open, SIGNAL(triggered()), this, SLOT(OnActionOpenTriggered()));
}

void ScatterPointGlyph::OnActionOpenTriggered() {
	if (scatter_point_data_ != NULL) {
		QMessageBox::warning(this, tr("Warning"), tr("Data has been opened!"));
		return;
	}

	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName("E:/Data/japan/Grid/timestep0_0003600.vtk");
	reader->SetReadAllScalars(1);
	reader->SetReadAllVectors(1);
	reader->Update();
	vtkDataObject* data = reader->GetOutput();
	scatter_point_data_ = vtkUnstructuredGrid::SafeDownCast(data);

	this->AddPointData2View();
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

	test_poly->SetPoints(scatter_point_data_->GetPoints());
	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
		vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(test_poly);
	vertexFilter->Update();
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->ShallowCopy(vertexFilter->GetOutput());

	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(1);
	colors->SetName("Colors");
	for (size_t i = 0; i < scatter_point_data_->GetPoints()->GetNumberOfPoints(); ++i){
		colors->InsertNextTuple1(128);
	}
	polydata->GetPointData()->SetScalars(colors);

	test_mapper->SetInputData(polydata);
	test_actor->SetMapper(test_mapper);
	test_actor->GetProperty()->SetPointSize(3);

	main_renderer_->AddActor(test_actor);
	main_renderer_->ResetCamera();
	main_view_->update();

	/// update the layer information
	RenderingLayer layer;
	layer.actor = test_actor;
	layer.name = "Point Layer";
	rendering_layer_vec_.push_back(layer);
}

void ScatterPointGlyph::contextMenuEvent(QContextMenuEvent *event){
	QCursor cur = this->cursor();
	context_menu_->exec(cur.pos());
}