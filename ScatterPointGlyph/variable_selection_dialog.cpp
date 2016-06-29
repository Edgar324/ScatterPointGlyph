#include "variable_selection_dialog.h"
#include <QtWidgets/QGridLayout>

VariableSelectionDialog::VariableSelectionDialog() {
	ui_.setupUi(this);

	connect(ui_.dim_red_checkbox, SIGNAL(stateChanged(int)), this, SLOT(OnAutoStateChanged(int)));
}

VariableSelectionDialog::~VariableSelectionDialog() {

}

void VariableSelectionDialog::SetDatasetInfo(int pnum, int vnum, std::vector<QString >& var_names) {
	ui_.point_num_label->setText(QString("%0").arg(pnum));
	ui_.variable_number_label->setText(QString("%0").arg(vnum));

	QGridLayout* grid_layout = new QGridLayout;
	for (int i = 0; i < var_names.size(); ++i) {
		QCheckBox* check_box = new QCheckBox(var_names[i]);
		check_box->setChecked(true);
		QDoubleSpinBox* spin_box = new QDoubleSpinBox;
		spin_box->setValue(1.0);
		spin_box->setRange(0, 1.0);

		var_check_boxes_.push_back(check_box);
		var_weight_spin_boxes_.push_back(spin_box);

		int row_index = i / 3;
		int colume_index = i % 3;
		grid_layout->addWidget(check_box, row_index, 2 * colume_index, 1, 1, Qt::AlignLeft);
		grid_layout->addWidget(spin_box, row_index, 2 * colume_index + 1, 1, 1, Qt::AlignLeft);
	}
	ui_.variable_selection_widget->setLayout(grid_layout);
}

void VariableSelectionDialog::GetSelectionResult(std::vector<float>& weights) {
	weights.resize(var_weight_spin_boxes_.size());
	for (int i = 0; i < var_check_boxes_.size(); ++i)
		if (var_check_boxes_[i]->isChecked())
			weights[i] = var_weight_spin_boxes_[i]->value();
		else
			weights[i] = -1.0;

	// normalize variable weight
	float accu_weight = 0;
	for (int i = 0; i < weights.size(); ++i)
		if (weights[i] >= 0) accu_weight += weights[i];
	if (accu_weight == 0) accu_weight = 1.0;
	for (int i = 0; i < weights.size(); ++i)
		if (weights[i] >= 0) weights[i] /= accu_weight;
}

void VariableSelectionDialog::OnAutoStateChanged(int is_enabled) {
	if (is_enabled == Qt::Checked)
		ui_.variable_selection_widget->setEnabled(false);
	else
		ui_.variable_selection_widget->setEnabled(true);
}