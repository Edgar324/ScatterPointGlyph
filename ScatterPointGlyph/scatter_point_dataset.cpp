#include "scatter_point_dataset.h"

ScatterPointDataset::ScatterPointDataset() {
	this->is_structured_data = false;
	this->w = 0;
	this->h = 0;
}

ScatterPointDataset::~ScatterPointDataset() {

}

void ScatterPointDataset::Sample(float left, float right, float bottom, float top) {
	point_pos.clear();
	point_values.clear();
	sample_index.clear();
	node_sample_map.clear();
	for (int i = 0; i < original_point_pos.size(); ++i)
		if ((original_point_pos[i][0] - left) * (original_point_pos[i][0] - right) <= 0  
			&& (original_point_pos[i][1] - bottom) * (original_point_pos[i][1] - top) <= 0) {
			point_pos.push_back(original_point_pos[i]);
			point_values.push_back(original_point_values[i]);
			sample_index.push_back(i);
			node_sample_map.insert(std::map< int, int >::value_type(i, (int)sample_index.size() - 1));
		}

	weights.resize(point_values[0].size());
	weights.assign(point_values[0].size(), 1.0 / point_values[0].size());

	NormalizePosition(this->point_pos, this->original_pos_ranges);
	NormalizeValues(this->point_values, this->original_value_ranges);
}

void ScatterPointDataset::DirectConstruct() {
	point_pos = original_point_pos;
	point_values = original_point_values;
	weights.resize(point_values[0].size());
	weights.assign(point_values[0].size(), 1.0 / point_values[0].size());
	sample_index.resize(point_pos.size());

	node_sample_map.clear();
	for (int i = 0; i < point_pos.size(); ++i) {
		sample_index[i] = i;
		node_sample_map.insert(std::map< int, int >::value_type(i, i));
	}

	NormalizePosition(this->point_pos, this->original_pos_ranges);
	NormalizeValues(this->point_values, this->original_value_ranges);
}

void ScatterPointDataset::NormalizeValues(std::vector< std::vector< float > >& vec, std::vector< std::vector< float > >& ranges){
	ranges.resize(vec[0].size());

	for (int i = 0; i < vec[0].size(); ++i){
		float minValue = 1e10;
		float maxValue = -1e10;

		for (int j = 0; j < vec.size(); ++j){
			if (minValue > vec[j][i]) minValue = vec[j][i];
			if (maxValue < vec[j][i]) maxValue = vec[j][i];
		}

		if (maxValue - minValue != 0) {
			for (int j = 0; j < vec.size(); ++j)
				vec[j][i] = (vec[j][i] - minValue) / (maxValue - minValue);
		}
		else {
			for (int j = 0; j < vec.size(); ++j) vec[j][i] = 0.5;
		}

		ranges[i].resize(2);
		ranges[i][0] = minValue;
		ranges[i][1] = maxValue;
	}
}

void ScatterPointDataset::NormalizePosition(std::vector< std::vector< float > >& vec, std::vector< std::vector< float > >& ranges) {
	ranges.resize(vec[0].size());

	float max_range = -1e10;
	for (int i = 0; i < vec[0].size(); ++i){
		float minValue = 1e10;
		float maxValue = -1e10;

		for (int j = 0; j < vec.size(); ++j){
			if (minValue > vec[j][i]) minValue = vec[j][i];
			if (maxValue < vec[j][i]) maxValue = vec[j][i];
		}
		if (maxValue - minValue > max_range) max_range = maxValue - minValue;

		ranges[i].resize(2);
		ranges[i][0] = minValue;
		ranges[i][1] = maxValue;
	}

	for (int i = 0; i < vec[0].size(); ++i){
		float mid_value = (ranges[i][0] + ranges[i][1]) / 2;
		
		ranges[i][0] = mid_value - max_range / 2;
		ranges[i][1] = mid_value + max_range / 2;

		for (int j = 0; j < vec.size(); ++j)
			vec[j][i] = (vec[j][i] - mid_value) / max_range + 0.5;
	}
}