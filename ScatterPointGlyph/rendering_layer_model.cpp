#include "rendering_layer_model.h"
#include <vtkActor.h>

RenderingLayerModel::RenderingLayerModel() 
	: max_id_(-1) {

}

RenderingLayerModel::~RenderingLayerModel() {

}

int RenderingLayerModel::AddLayer(QString& name, vtkActor* layer_actor, bool is_visible) {
	max_id_++;

	RenderingLayer* layer = new RenderingLayer(max_id_);
	layer->name = name;
	layer->layer_actor = layer_actor;
	layer->layer_actor->SetVisibility(is_visible);

	rendering_layer_list_.push_back(layer);

	emit ModelChanged();

	return max_id_;
}

bool RenderingLayerModel::RemoveLayer(int id) {
	std::list< RenderingLayer* >::iterator iter = rendering_layer_list_.begin();
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
	std::list< RenderingLayer* >::iterator iter = rendering_layer_list_.begin();
	while (iter != rendering_layer_list_.end()) {
		layers.push_back(*iter);
		iter++;
	}
}

RenderingLayer* RenderingLayerModel::GetLayer(int index) {
	if (index < this->rendering_layer_list_.size()) {
		std::list< RenderingLayer* >::iterator iter = this->rendering_layer_list_.begin();
		int temp_index = 0;
		while (temp_index < index && iter != this->rendering_layer_list_.end()) {
			temp_index++;
			iter++;
		}
		return (*iter);
	} else {
		return NULL;
	}
}

void RenderingLayerModel::SetLayerVisibility(int index, bool visibility) {
	RenderingLayer* layer = this->GetLayer(index);
	if (layer != NULL) {
		layer->layer_actor->SetVisibility(visibility);
		emit LayerPropertyChanged();
	}
}

void RenderingLayerModel::SetLayerName(int index, QString name) {
	RenderingLayer* layer = this->GetLayer(index);
	if (layer != NULL) {
		layer->name = name;
	}
}