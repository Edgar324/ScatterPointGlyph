#ifndef VARIABLE_SELECTION_DIALOG_H_
#define VARIABLE_SELECTION_DIALOG_H_

#include <QtWidgets/QDialog>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QCheckBox>
#include <vector>
#include "ui_variable_selection_dialog.h"

class VariableSelectionDialog : public QDialog
{
	Q_OBJECT

public:
	VariableSelectionDialog();
	~VariableSelectionDialog();

	void SetDatasetInfo(int pnum, int vnum, std::vector<QString>& var_names);
	void GetSelectionResult(std::vector<float>& weights);
	bool IsAutomaticDimReduction() { return ui_.dim_red_checkbox->isChecked();  }
	int GetDimNumber() { return ui_.auto_dim_number->value(); }

private:
	Ui::VarSelectionDialog ui_;
	std::vector<QCheckBox*> var_check_boxes_;
	std::vector<QDoubleSpinBox*> var_weight_spin_boxes_;

	private slots:
	void OnAutoStateChanged(int);
};

#endif