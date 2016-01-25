#include "layer_control_widget.h"

#include <QtWidgets/QTableView>
#include <QtGui/QStandardItemModel>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QGroupBox>
#include <vtk3DWidget.h>
#include "rendering_layer_model.h"

LayerControlWidget::LayerControlWidget()
	: layer_model_(NULL) {
	this->InitWidget();
}

LayerControlWidget::~LayerControlWidget() {

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

	QGroupBox* layer_groupbox = new QGroupBox(QString("Rendering Layers"));
	QVBoxLayout* layout2 = new QVBoxLayout;
	layout2->addWidget(rendering_layer_tableview_);
	layer_groupbox->setLayout(layout2);

	QVBoxLayout* main_layout = new QVBoxLayout;
	main_layout->addWidget(layer_groupbox);
	this->setLayout(main_layout);

	connect(layer_item_model_, SIGNAL(itemChanged(QStandardItem*)), this, SLOT(OnItemModelChanged(QStandardItem*)));
}

void LayerControlWidget::UpdateWidget() {
	if (layer_model_ != NULL) {
		disconnect(layer_item_model_, SIGNAL(itemChanged(QStandardItem*)), this, SLOT(OnItemModelChanged(QStandardItem*)));

		std::vector< RenderingLayer* > layers;
		layer_model_->GetAllLayers(layers);
		layer_item_model_->setRowCount(layers.size());
		for (size_t i = 0; i < layers.size(); ++i) {
			layer_item_model_->setData(layer_item_model_->index(i, 0), layers[i]->layer_id(), Qt::DisplayRole);
			layer_item_model_->setData(layer_item_model_->index(i, 1), layers[i]->layer_widget->GetEnabled() ? Qt::Checked : Qt::Unchecked, Qt::CheckStateRole);
			layer_item_model_->setData(layer_item_model_->index(i, 2), layers[i]->name, Qt::DisplayRole);
			layer_item_model_->item(i, 1)->setCheckable(true);
			layer_item_model_->item(i, 2)->setTextAlignment(Qt::AlignCenter);
		}
		rendering_layer_tableview_->resizeColumnsToContents();

		connect(layer_item_model_, SIGNAL(itemChanged(QStandardItem*)), this, SLOT(OnItemModelChanged(QStandardItem*)));
	}
}

void LayerControlWidget::OnLayerModelChanged() {
	this->UpdateWidget();
}

void LayerControlWidget::OnItemModelChanged(QStandardItem* item) {
	int row_index = item->row();
	int colum_index = item->column();

	if (colum_index == 1) {
		layer_model_->SetLayerVisibility(row_index, item->checkState());
	} else if (colum_index == 2) {
		layer_model_->SetLayerName(row_index, item->text());
	}
}