#include "layer_control_widget.h"

#include <QtWidgets/QTableView>
#include <QtGui/QStandardItemModel>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QGroupBox>
#include <vtkActor.h>
#include "rendering_layer_model.h"

LayerControlWidget::LayerControlWidget()
	: layer_model_(NULL) {
	this->InitWidget();
}

LayerControlWidget::~LayerControlWidget() {

}

void LayerControlWidget::SetDatasetInfo(QString& name, int point_num) {
	dataset_name_ = name;
	dataset_point_number_ = point_num;
}

void LayerControlWidget::SetRenderingLayerModel(RenderingLayerModel* model) {
	if (layer_model_ != NULL) {
		delete layer_model_;
		layer_model_ = NULL;
	}
	layer_model_ = model;
	if (layer_model_ != NULL) {
		connect(layer_model_, SIGNAL(ModelChanged()), this, SLOT(OnLayerModelChanged()));
	}
}

void LayerControlWidget::InitWidget() {
	dataset_tableview_ = new QTableView;
	dataset_item_model_ = new QStandardItemModel;
	QStringList dataset_headers;
	dataset_headers << "Name" << "Point Count";
	dataset_item_model_->setHorizontalHeaderLabels(dataset_headers);
	dataset_tableview_->setModel(dataset_item_model_);
	dataset_tableview_->resizeColumnsToContents();

	rendering_layer_tableview_ = new QTableView;
	layer_item_model_ = new QStandardItemModel;
	QStringList layer_headers;
	layer_headers << "ID" << "Visibility" << "Name";
	layer_item_model_->setHorizontalHeaderLabels(layer_headers);
	rendering_layer_tableview_->setModel(layer_item_model_);
	rendering_layer_tableview_->setColumnHidden(0, true);
	DisableColumnEidtDelegate *dced = new DisableColumnEidtDelegate();
	rendering_layer_tableview_->setItemDelegateForColumn(1, dced);
	rendering_layer_tableview_->resizeColumnsToContents();

	
	QGroupBox* dataset_groupbox = new QGroupBox(QString("Dataset"));
	QVBoxLayout* layout1 = new QVBoxLayout;
	layout1->addWidget(dataset_tableview_);
	dataset_groupbox->setLayout(layout1);

	QGroupBox* layer_groupbox = new QGroupBox(QString("Rendering Layers"));
	QVBoxLayout* layout2 = new QVBoxLayout;
	layout2->addWidget(rendering_layer_tableview_);
	layer_groupbox->setLayout(layout2);

	QVBoxLayout* main_layout = new QVBoxLayout;
	main_layout->addWidget(dataset_groupbox);
	main_layout->addWidget(layer_groupbox);
	this->setLayout(main_layout);
}

void LayerControlWidget::UpdateWidget() {
	if (dataset_name_.length() != 0) {
		dataset_item_model_->setRowCount(1);
		dataset_item_model_->setData(dataset_item_model_->index(0, 0), dataset_name_, Qt::DisplayRole);
		dataset_item_model_->setData(dataset_item_model_->index(0, 1), dataset_point_number_, Qt::DisplayRole);
	}
	else {
		dataset_item_model_->setRowCount(0);
	}

	if (layer_model_ != NULL) {
		std::vector< RenderingLayer* > layers;
		layer_model_->GetAllLayers(layers);
		layer_item_model_->setRowCount(layers.size());
		for (size_t i = 0; i < layers.size(); ++i) {
			layer_item_model_->setData(layer_item_model_->index(i, 0), layers[i]->layer_id(), Qt::DisplayRole);
			layer_item_model_->setData(layer_item_model_->index(i, 1), layers[i]->layer_actor->GetVisibility() ? Qt::Checked : Qt::Unchecked, Qt::CheckStateRole);
			layer_item_model_->setData(layer_item_model_->index(i, 2), layers[i]->name, Qt::DisplayRole);
		}
	}
}

void LayerControlWidget::OnLayerModelChanged() {
	this->UpdateWidget();
}