#include "scatter_point_dataset.h"

ScatterPointDataset::ScatterPointDataset() {

}

ScatterPointDataset::~ScatterPointDataset() {

}

void ScatterPointDataset::Sample(int point_num, float left, float right, float bottom, float top) {

}

void ScatterPointDataset::DirectConstruct() {
	point_pos = original_point_pos;

	point_values.resize(original_point_values.size());
	for (int i = 0; i < original_point_pos.size(); ++i) {
		float temp_value = 0;
		for (int j = 0; j < weights.size(); ++j)
			temp_value += original_point_values[i][j] * weights[j];
		point_values[i] = temp_value;
	}
}