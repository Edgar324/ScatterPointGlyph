#include "hier_para_widget.h"

HierParaWidget::HierParaWidget() {
	ui_.setupUi(this);

	this->SetMaxClusterNumber(1);

	connect(ui_.cluster_number_edit, SIGNAL(editingFinished()), this, SLOT(OnEditValueChanged()));
	connect(ui_.cluster_number_slider, SIGNAL(sliderReleased()), this, SLOT(OnSliderValueChanged()));
}

HierParaWidget::~HierParaWidget() {

}

void HierParaWidget::SetMaxClusterNumber(int num) {
	max_cluster_num_ = num;

	ui_.cluster_number_slider->setRange(1, max_cluster_num_);
	ui_.cluster_number_slider->setSingleStep(1);
	ui_.cluster_number_slider->setValue(max_cluster_num_);

	ui_.cluster_number_edit->setText(QString::number(max_cluster_num_));
}

void HierParaWidget::OnSliderValueChanged() {
	ui_.cluster_number_edit->setText(QString::number(ui_.cluster_number_slider->value()));

	emit ClusterNumberChanged(ui_.cluster_number_slider->value());
}

void HierParaWidget::OnEditValueChanged() {
	ui_.cluster_number_slider->setValue(ui_.cluster_number_edit->text().toInt());

	emit ClusterNumberChanged(ui_.cluster_number_edit->text().toInt());
}