#include "rendering_layer_model.h"

RenderingLayerModel::RenderingLayerModel() 
	: max_id_(-1) {

}

RenderingLayerModel::~RenderingLayerModel() {

}

int RenderingLayerModel::AddLayer(QString& name, vtkActor* layer_actor) {
	max_id_++;

	RenderingLayer* layer = new RenderingLayer(max_id_);
	layer->name = name;
	layer->layer_actor = layer_actor;

	rendering_layer_list_.push_back(layer);

	emit ModelChanged();

	return max_id_;
}

bool RenderingLayerModel::RemoveLayer(int id) {
	std::list< RenderingLayer* >::iterator iter;
	while (iter != rendering_layer_list_.end() && (*iter)->layer_id() != id) iter++;
	if (iter != rendering_layer_list_.end()) {
		rendering_layer_list_.erase(iter);

		emit ModelChanged();

		return true;
	} else 
		return false;
}

void RenderingLayerModel::GetAllLayers(std::vector< RenderingLayer* >& layers) {
	layers.clear();
	std::list< RenderingLayer* >::iterator iter;
	while (iter != rendering_layer_list_.end()) {
		layers.push_back(*iter);
		iter++;
	}
}