#ifndef HIER_PARA_WIDGET_H_
#define HIER_PARA_WIDGET_H_

#include <QtWidgets/QWidget>
#include "ui_hier_para_widget.h"

class HierParaWidget : public QWidget
{
	Q_OBJECT

public:
	HierParaWidget();
	~HierParaWidget();

	void SetMaxClusterNumber(int num);

signals:
	void ClusterNumberChanged(int);

private:
	Ui::HierParaWidget ui_;

	int max_cluster_num_;

private slots:
	void OnSliderValueChanged();
	void OnEditValueChanged();
};

#endif