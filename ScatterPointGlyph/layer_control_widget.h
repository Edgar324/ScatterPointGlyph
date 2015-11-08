#ifndef LAYER_CONTROL_WIDGET_H_
#define LAYER_CONTROL_WIDGET_H_

#include <QtWidgets/QWidget>

class RenderingLayerModel;
class QTableView;
class QStandardItemModel;

class LayerControlWidget : public QWidget
{
	Q_OBJECT

public:
	LayerControlWidget();
	~LayerControlWidget();

	void SetDatasetInfo(QString& name, int point_num);
	void SetRenderingLayerModel(RenderingLayerModel* model);

private:
	RenderingLayerModel* layer_model_;

	QString dataset_name_;
	int dataset_point_number_;

	QTableView* dataset_tableview_;
	QStandardItemModel* dataset_item_model_;
	QTableView* rendering_layer_tableview_;
	QStandardItemModel* layer_item_model_;

	void InitWidget();
	void UpdateWidget();

private slots:
	void OnLayerModelChanged();
};

#endif