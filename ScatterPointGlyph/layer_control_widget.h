#ifndef LAYER_CONTROL_WIDGET_H_
#define LAYER_CONTROL_WIDGET_H_

#include <QtWidgets/QWidget>
#include <QStyledItemDelegate>

class RenderingLayerModel;
class QTableView;
class QStandardItem;
class QStandardItemModel;

class DisableColumnEidtDelegate : public QStyledItemDelegate
{
	Q_OBJECT
public:
	DisableColumnEidtDelegate(QObject *parent = 0) : QStyledItemDelegate(parent) { }
	QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
		const QModelIndex &index) const {
		return NULL;
	}
};


class LayerControlWidget : public QWidget
{
	Q_OBJECT

public:
	LayerControlWidget();
	~LayerControlWidget();

	void SetRenderingLayerModel(RenderingLayerModel* model);

private:
	RenderingLayerModel* layer_model_;

	QTableView* rendering_layer_tableview_;
	QStandardItemModel* layer_item_model_;

	void InitWidget();
	void UpdateWidget();

private slots:
	void OnLayerModelChanged();
	void OnItemModelChanged(QStandardItem* item);
};

#endif