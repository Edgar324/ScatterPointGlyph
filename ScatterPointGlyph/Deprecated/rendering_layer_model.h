#ifndef RENDERING_LAYER_MODEL_H_
#define RENDERING_LAYER_MODEL_H_

#include <list>
#include <string>
#include <vector>
#include <QtCore/QString>
#include <QtCore/QObject>

class vtk3DWidget;

class RenderingLayer
{
public:
	RenderingLayer(int id) {
		layer_id_ = id;
		name = "";
		layer_widget = NULL;
	}
	~RenderingLayer() {}

	QString name;
	vtk3DWidget* layer_widget;

	int layer_id() { return layer_id_; }

private:
	int layer_id_;
};

class RenderingLayerModel : public QObject
{
	Q_OBJECT

public:
	RenderingLayerModel();
	~RenderingLayerModel();

	int AddLayer(QString& name, vtk3DWidget* layer_actor, bool is_visible = true);
	bool RemoveLayer(int id);
	
	void GetAllLayers(std::vector< RenderingLayer* >& layers);
	RenderingLayer* GetLayer(int index);

	void SetLayerVisibility(int index, bool visibility);
	void SetLayerName(int index, QString name);

signals:
	void ModelChanged();
	void LayerPropertyChanged();

private:
	int max_id_;
	std::list< RenderingLayer* > rendering_layer_list_;
};

#endif