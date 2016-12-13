/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "image_widget.h"
#include <vtkCellPicker.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>

vtkStandardNewMacro(ImageWidget);

ImageWidget::ImageWidget() {
    this->EventCallbackCommand->SetCallback(ImageWidget::ProcessEvents);
    this->picker_ = vtkCellPicker::New();
    this->picker_->PickFromListOn();
}

ImageWidget::~ImageWidget() {

}

void ImageWidget::SetData(int w, int h, vector<float>& rgba_values) {
    image_width_ = w;
    image_height_ = h;
    image_rgba_values_ = rgba_values;

    this->BuildRepresentation();
}

void ImageWidget::SetEnabled(int enabling) {
    if (!this->Interactor) {
		vtkErrorMacro(<< "The interactor must be set prior to enabling/disabling widget");
		return;
	}

    if (this->image_actor_ == NULL) {
        vtkErrorMacro(<< "The image must be set prior to enabling/disabling widget");
		return;
    }

	if (enabling) {
		if (this->Enabled) return;

        if (!this->CurrentRenderer) {
            this->SetCurrentRenderer(this->Interactor->FindPokedRenderer(
                this->Interactor->GetLastEventPosition()[0],
                this->Interactor->GetLastEventPosition()[1]));
            if (this->CurrentRenderer == NULL) return;
        }

		this->Enabled = 1;

        this->Interactor->AddObserver(vtkCommand::LeftButtonPressEvent,
            this->EventCallbackCommand, this->Priority);
        this->Interactor->AddObserver(vtkCommand::MouseMoveEvent,
            this->EventCallbackCommand, this->Priority);
        this->Interactor->AddObserver(vtkCommand::LeftButtonReleaseEvent,
            this->EventCallbackCommand, this->Priority);
        this->Interactor->AddObserver(vtkCommand::RightButtonPressEvent,
            this->EventCallbackCommand, this->Priority);
        this->Interactor->AddObserver(vtkCommand::RightButtonReleaseEvent,
            this->EventCallbackCommand, this->Priority);

		this->CurrentRenderer->AddActor(this->image_actor_);

        this->InvokeEvent(vtkCommand::EnableEvent,NULL);
	}
	else {
		if (!this->Enabled) return;

		this->Enabled = 0;

        this->Interactor->RemoveObserver(this->EventCallbackCommand);

		this->CurrentRenderer->RemoveActor(this->image_actor_);

        this->InvokeEvent(vtkCommand::DisableEvent,NULL);
        this->SetCurrentRenderer(NULL);
	}
}

void ImageWidget::Modified() {

}

void ImageWidget::BuildRepresentation() {
    vtkSmartPointer<vtkImageData> image_data = vtkSmartPointer<vtkImageData>::New();
    image_data->SetDimensions(image_width_, image_height_, 1);
    image_data->AllocateScalars(VTK_FLOAT, 4);

    int accu_index = 0;
    for (int y = 0; y < image_height_; y++)
        for (int x = 0; x < image_width_; x++) {
            float* pixel = static_cast<float*>(image_data->GetScalarPointer(x, y, 0));
            pixel[0] = image_rgba_values_[accu_index];
            pixel[1] = image_rgba_values_[accu_index + 1];
            pixel[2] = image_rgba_values_[accu_index + 2];
            pixel[3] = image_rgba_values_[accu_index + 3];
        }

    if (this->image_actor_ == NULL) {
        this->image_actor_ = vtkImageActor::New();
    }

    this->image_actor_->SetInputData(image_data);
    this->image_actor_->Modified();
}

void ImageWidget::ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata) {
    ImageWidget* self = reinterpret_cast<ImageWidget*>(clientdata);

	switch (event)
	{
	case vtkCommand::LeftButtonPressEvent:
		self->OnLeftButtonDown();
		break;
	case vtkCommand::LeftButtonReleaseEvent:
		self->OnLeftButtonUp();
		break;
	case vtkCommand::RightButtonPressEvent:
		self->OnRightButtonDown();
		break;
	case vtkCommand::RightButtonReleaseEvent:
		self->OnRightButtonUp();
		break;
	case vtkCommand::MouseMoveEvent:
		self->OnMouseMove();
		break;
	default:
		break;
	}
}

void ImageWidget::OnMouseMove() {
    cout << "Image widget move" << endl;
}

void ImageWidget::OnLeftButtonDown() {

}

void ImageWidget::OnLeftButtonUp() {

}

void ImageWidget::OnRightButtonDown() {

}

void ImageWidget::OnRightButtonUp() {

}