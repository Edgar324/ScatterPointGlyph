#include "scatter_point_dataset.h"
#include <assert.h>
#include "color_mapping_generator.h"
#include "mds_projector.h"
#include "pca_projector.h"
#include "tsne_projector.h"

ScatterPointDataset::ScatterPointDataset() {
	
}

ScatterPointDataset::~ScatterPointDataset() {

}

void ScatterPointDataset::ClearData() {
	var_num = 0;
	point_num = 0;
	normalized_point_pos.clear();
	normalized_point_values.clear();

	original_point_pos.clear();
	original_pos_ranges.clear();

	var_names.clear();
	var_weights.clear();
    var_colors.clear();
	original_point_values.clear();
	original_value_ranges.clear();
}

void ScatterPointDataset::DirectConstruct() {
	var_num = var_names.size();
	point_num = original_point_values[0].size();

    is_valid.resize(point_num, true);

    var_weights.resize(var_num, 1.0 / var_num);

    // assign variable colors
    ColorMappingGenerator::GetInstance()->GetQualitativeColors(var_num, var_colors);

	normalized_point_pos = original_point_pos;
	normalized_point_values = original_point_values;

	NormalizePosition(this->normalized_point_pos, this->original_pos_ranges);
	NormalizeValues(this->normalized_point_values, this->original_value_ranges);
}

void ScatterPointDataset::ConstructWithValueRanges(vector<vector<float>>& ranges) {
    var_num = var_names.size();
	point_num = original_point_values[0].size();
    is_valid.resize(point_num, true);

    var_weights.resize(var_num, 1.0 / var_num);

    // assign variable colors
    ColorMappingGenerator::GetInstance()->GetQualitativeColors(var_num, var_colors);

	normalized_point_values = original_point_values;
    this->original_value_ranges = ranges;
    for (int i = 0; i < point_num; ++i) {
        is_valid[i] = true;
        for (int j = 0; j < var_num; ++j) 
            is_valid[i] = is_valid[i] && ((original_point_values[j][i] - ranges[j][0]) * (original_point_values[j][i] - ranges[j][1]) <= 0);
        /*if (is_valid[i]) {
            for (int j = 0; j < var_num; ++j)
                normalized_point_values[j][i] = (original_point_values[j][i] - ranges[j][0]) / (ranges[j][1] - ranges[j][0]);
        }*/
    }

    normalized_point_pos = original_point_pos;
	normalized_point_values = original_point_values;

	NormalizePosition(this->normalized_point_pos, this->original_pos_ranges);
	NormalizeValues(this->normalized_point_values, this->original_value_ranges);
}

void ScatterPointDataset::NormalizePos() {
    normalized_point_pos = original_point_pos;

	NormalizePosition(this->normalized_point_pos, this->original_pos_ranges);
}

void ScatterPointDataset::ApplyMds() {
	MdsProjector projector;
    projector.Project(this->original_point_values, this->original_point_pos);

    this->DirectConstruct();
}

void ScatterPointDataset::ApplyTsne() {
    TsneProjector projector;
    projector.Project(this->original_point_values, this->original_point_pos);

    this->DirectConstruct();
}

void ScatterPointDataset::ApplyNormal(int axis_one, int axis_two) {
    int point_num = this->original_point_values[0].size();
    this->original_point_pos.resize(2);
    this->original_point_pos[0].resize(point_num);
    this->original_point_pos[1].resize(point_num);

    is_valid.resize(point_num, true);
    normalized_point_values = original_point_values;
    NormalizeValues(this->normalized_point_values, this->original_value_ranges);

    for (int i = 0; i < point_num; i++) {
        this->original_point_pos[0][i] = this->normalized_point_values[axis_one][i] * 1000;
        this->original_point_pos[1][i] = this->normalized_point_values[axis_two][i] * 1000;
    }

    this->DirectConstruct();
}

void ScatterPointDataset::NormalizeValues(vector<vector<float>>& vec, vector<vector<float>>& ranges){
	ranges.resize(vec.size());

	for (int i = 0; i < vec.size(); ++i){
		float minValue = 1e10;
		float maxValue = -1e10;

		for (int j = 0; j < vec[i].size(); ++j){
            if (!is_valid[j]) continue;
			if (minValue > vec[i][j]) minValue = vec[i][j];
			if (maxValue < vec[i][j]) maxValue = vec[i][j];
		}

		if (maxValue - minValue != 0) {
            for (int j = 0; j < vec[i].size(); ++j) {
                if (!is_valid[j]) continue;
                vec[i][j] = (vec[i][j] - minValue) / (maxValue - minValue);
            }
				
		}
		else {
			for (int j = 0; j < vec[i].size(); ++j) vec[i][j] = 0.5;
		}

		ranges[i].resize(2);
		ranges[i][0] = minValue;
		ranges[i][1] = maxValue;
	}
}

void ScatterPointDataset::NormalizePosition(vector<vector<float>>& vec, vector<vector<float>>& ranges) {
	ranges.resize(vec.size());

	max_pos_range = -1e10;
	for (int i = 0; i < vec.size(); ++i){
		float minValue = 1e10;
		float maxValue = -1e10;

		for (int j = 0; j < vec[i].size(); ++j){
            if (!is_valid[j]) continue;
			if (minValue > vec[i][j]) minValue = vec[i][j];
			if (maxValue < vec[i][j]) maxValue = vec[i][j];
		}
		if (maxValue - minValue > max_pos_range) max_pos_range = maxValue - minValue;

		ranges[i].resize(2);
		ranges[i][0] = minValue;
		ranges[i][1] = maxValue;
	}

	for (int i = 0; i < vec.size(); ++i){
		float mid_value = (ranges[i][0] + ranges[i][1]) / 2;
		
		ranges[i][0] = mid_value - max_pos_range / 2;
		ranges[i][1] = mid_value + max_pos_range / 2;

        for (int j = 0; j < vec[i].size(); ++j) {
            if (!is_valid[j]) continue;
			vec[i][j] = (vec[i][j] - mid_value) / max_pos_range + 0.5;
        }
	}
}
