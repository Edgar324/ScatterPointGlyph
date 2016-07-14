#include "qvtk_rendering_widget.h"
#include <iostream>
#include <QTimer>
#include <QWheelEvent>
#include <QMouseEvent>
#include <QDragMoveEvent>
#include <QtWidgets/QToolTip>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleImage.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include "cnode.h"

QVtkRenderingWidget::QVtkRenderingWidget() 
	: QVTKWidget(), is_wheel_updated_(false) {
	timer_ = new QTimer;
	timer_->setSingleShot(true);

	this->setFocusPolicy(Qt::StrongFocus);

    main_renderer_ = vtkRenderer::New();
	this->GetRenderWindow()->AddRenderer(main_renderer_);
	main_renderer_->SetViewport(0.0, 0.0, 1.0, 1.0);
	main_renderer_->SetBackground(1.0, 1.0, 1.0);

    indicator_renderer_ = vtkRenderer::New();
    indicator_renderer_->SetViewport(0.85, 0.8, 1.0, 1.0);
	indicator_renderer_->SetBackground(1.0, 1.0, 1.0);
    this->GetRenderWindow()->AddRenderer(indicator_renderer_);

    main_renderer_->GetActiveCamera()->SetParallelProjection(0);

	vtkSmartPointer<vtkInteractorStyleImage> imageStyle =
		vtkSmartPointer<vtkInteractorStyleImage>::New();
	this->GetRenderWindow()->GetInteractor()->SetInteractorStyle(imageStyle);
	this->GetRenderWindow()->SetPointSmoothing(1);
    this->GetRenderWindow()->SetLineSmoothing(1);

    this->selection_brush_poly = vtkPolyData::New();
	this->selection_brush_mapper = vtkPolyDataMapper::New();
	this->selection_brush_actor = vtkActor::New();
	this->selection_brush_mapper->SetInputData(this->selection_brush_poly);
	this->selection_brush_actor->SetMapper(this->selection_brush_mapper);
	this->selection_brush_actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
	this->selection_brush_actor->GetProperty()->SetLineWidth(3.0);

    this->main_renderer_->AddActor(selection_brush_actor);

	connect(timer_, SIGNAL(timeout()), this, SLOT(OnTimeout()));
}

QVtkRenderingWidget::~QVtkRenderingWidget() {

}

void QVtkRenderingWidget::SetBrushEnabled(bool enabled) {
    if (!enabled) {
        this->selection_brush_poly->Initialize();
		vtkPoints* points = vtkPoints::New();
		this->selection_brush_poly->SetPoints(points);
		vtkCellArray* line_array = vtkCellArray::New();
		this->selection_brush_poly->SetLines(line_array);
		this->selection_brush_poly->Modified();
		this->selection_brush_actor->Modified();

        brushing_path.clear();

        this->update();
    }

    this->is_brush_enabled_ = enabled;
}

float QVtkRenderingWidget::GetDistancePerPixel() {
    double point_one[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 1, 0, 0, point_one);
	double point_two[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 2, 0, 0, point_two);

	return (1.0 / abs(point_one[0] - point_two[0]));
}

void QVtkRenderingWidget::GetViewport(float& left, float& right, float& bottom, float& top) {
    double viewport[4];
	this->main_renderer_->GetViewport(viewport);
	double point_one[4];
	vtkInteractorObserver::ComputeDisplayToWorld(this->main_renderer_, viewport[0] * this->width(), viewport[1] * this->height(), 0, point_one);
	double point_two[4];
	vtkInteractorObserver::ComputeDisplayToWorld(this->main_renderer_, viewport[2] * this->width(), viewport[3] * this->height(), 0, point_two);
	left = point_one[0];
	right = point_two[0];
	bottom = point_one[1];
	top = point_two[1];
}

void QVtkRenderingWidget::wheelEvent(QWheelEvent* event) {
	is_wheel_updated_ = true;
	if (!timer_->isActive()) {
		timer_->start(200);
	}

	QVTKWidget::wheelEvent(event);
}

void QVtkRenderingWidget::mousePressEvent(QMouseEvent* event) {
    if (is_brush_enabled_) {
        this->selection_brush_poly->Initialize();
		vtkPoints* points = vtkPoints::New();
		this->selection_brush_poly->SetPoints(points);
		vtkCellArray* line_array = vtkCellArray::New();
		this->selection_brush_poly->SetLines(line_array);
		this->selection_brush_poly->Modified();
		this->selection_brush_actor->Modified();

        brushing_path.clear();

        this->update();
    } else {
        QVTKWidget::mousePressEvent(event);
    }
}

void QVtkRenderingWidget::mouseMoveEvent(QMouseEvent* event) {
	global_mouse_pos_ = event->globalPos();

    if (is_brush_enabled_ && event->buttons() & Qt::LeftButton) {
        double display_pos[4];
		display_pos[0] = event->x();
		display_pos[1] = this->height() - event->y();
		display_pos[2] = 0;
		display_pos[3] = 1;
		this->main_renderer_->SetDisplayPoint(display_pos);
		this->main_renderer_->DisplayToWorld();
		double* world_pos = this->main_renderer_->GetWorldPoint();
        brushing_path.push_back(world_pos[0]);
        brushing_path.push_back(world_pos[1]);

        this->selection_brush_poly->Initialize();
		vtkPoints* points = vtkPoints::New();
		this->selection_brush_poly->SetPoints(points);
		vtkCellArray* line_array = vtkCellArray::New();
		this->selection_brush_poly->SetLines(line_array);

        for (int i = 0; i < brushing_path.size() / 2; ++i) {
            points->InsertNextPoint(brushing_path[2 * i], brushing_path[2 * i + 1], 0.001);
        }

		int pnum = this->selection_brush_poly->GetPoints()->GetNumberOfPoints();
		vtkIdType cell_ids[2];
		for (int i = 0; i < pnum; ++i){
			cell_ids[0] = i;
			cell_ids[1] = (i + 1) % pnum;
			this->selection_brush_poly->InsertNextCell(VTK_LINE, 2, cell_ids);
		}

		this->selection_brush_poly->Modified();
		this->selection_brush_actor->Modified();

        this->update();
    } else {
        QVTKWidget::mouseMoveEvent(event);
    }
}

void QVtkRenderingWidget::mouseReleaseEvent(QMouseEvent* event) {
    if (is_brush_enabled_) {
        emit SelectionChanged();

        this->selection_brush_poly->Initialize();
		vtkPoints* points = vtkPoints::New();
		this->selection_brush_poly->SetPoints(points);
		vtkCellArray* line_array = vtkCellArray::New();
		this->selection_brush_poly->SetLines(line_array);
		this->selection_brush_poly->Modified();
		this->selection_brush_actor->Modified();

        brushing_path.clear();

        this->update();
    } else {
        QVTKWidget::mouseReleaseEvent(event);
    }
}

void QVtkRenderingWidget::OnTimeout() {
	if (is_wheel_updated_) {
		is_wheel_updated_ = false;
		timer_->start(200);
	} else {
		emit ViewUpdated();
	}
}

void QVtkRenderingWidget::ShowTooltip(QString tip) {
	/*QString tip_str = "";
	tip_str = "Cluster Information: <br>";
	tip_str += QString("Point Count: %0<br>").arg(point_count);
	tip_str += "Variable Values: Mean(SD)<br>";
	QString temp_str = "";
	temp_str += "<font color=red>";
	temp_str += axis_name;
	temp_str += QString(":  %0(%1)").arg(average_value).arg(variance_value);
	temp_str += "</font><br>";
	tip_str += temp_str;*/
	QToolTip::showText(global_mouse_pos_, tip);
}

void QVtkRenderingWidget::HideTooltip() {
	QToolTip::hideText();
}
