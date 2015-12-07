#include "scatter_point_dataset.h"

ScatterPointDataset::ScatterPointDataset() {
	this->is_structured_data = false;
}

ScatterPointDataset::~ScatterPointDataset() {

}

void ScatterPointDataset::Sample(int point_num, float left, float right, float bottom, float top) {

}

void ScatterPointDataset::DirectConstruct() {
	point_pos = original_point_pos;
	point_values = original_point_values;
}