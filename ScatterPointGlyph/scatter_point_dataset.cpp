#include "scatter_point_dataset.h"
#include <assert.h>
#include "SimpleMatrix.h"
#include "color_mapping_generator.h"

ScatterPointDataset::ScatterPointDataset() {
	
}

ScatterPointDataset::~ScatterPointDataset() {

}

void ScatterPointDataset::ClearData() {
	var_num = 0;
	point_num = 0;
	point_pos.clear();
	point_values.clear();

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
	point_num = original_point_values.size();

    // assign variable colors
    ColorMappingGenerator::GetInstance()->GetQualitativeColors(var_num, var_colors);

	point_pos = original_point_pos;
	point_values = original_point_values;

	NormalizePosition(this->point_pos, this->original_pos_ranges);
	NormalizeValues(this->point_values, this->original_value_ranges);

	// construct adaptive rate
	std::vector< std::vector< float > > point_dis;
	point_dis.resize(this->point_pos.size());
	for (int i = 0; i < point_num; ++i)
		point_dis[i].resize(point_pos.size());
	for (int i = 0; i < this->point_pos.size() - 1; ++i) {
		point_dis[i][i] = 0;
		for (int j = i + 1; j < this->point_pos.size(); ++j) {
			float dis = sqrt(pow(point_pos[i][0] - point_pos[j][0], 2) + pow(point_pos[i][1] - point_pos[j][1], 2));
			point_dis[i][j] = dis;
			point_dis[j][i] = dis;
		}
	}
	for (int i = 0; i < point_dis.size(); ++i)
		sort(point_dis[i].begin(), point_dis[i].end());
	adaptive_rate.resize(point_num);
	int kNum = 20;
	for (int i = 0; i < point_num; ++i) {
		float ave_dis = 0;
		for (int j = 0; j < kNum; ++j)
			ave_dis += point_dis[i][j];
		ave_dis /= (kNum - 1);
		if (ave_dis > 0.1) ave_dis = 0.1;
		adaptive_rate[i] = ave_dis;
	}
}

void ScatterPointDataset::AutoDimReduction(int dim_num) {

}

void ScatterPointDataset::ManualSelectDim(std::vector< bool >& is_dim_selected) {
    selected_vars.clear();
    for (int i = 0; i < is_dim_selected.size(); ++i)
        if (is_dim_selected[i]) selected_vars.push_back(i);
        else var_weights[i] = 0.0;

	/*int accu_num = 0;
	for (int i = 0; i < is_dim_selected.size(); ++i) {
		if (is_dim_selected[i]) {
			for (int j = 0; j < original_point_values.size(); ++j)
				original_point_values[j][accu_num] = original_point_values[j][i];
			var_names[accu_num] = var_names[i];
			accu_num++;
		}
	}

	for (int i = 0; i < original_point_values.size(); ++i)
		original_point_values[i].resize(accu_num);
	var_names.resize(accu_num);*/
}

void ScatterPointDataset::ExecMds() {
	smat::Matrix<double> *X0 = new smat::Matrix<double>(original_point_values.size(), selected_vars.size(), 0.0);

	for (int i = 0; i < selected_vars.size(); ++i) {
		float min_value = 1e30;
		float max_value = -1e30;
		for (int j = 0; j < original_point_values.size(); ++j) {
			if (original_point_values[j][selected_vars[i]] > max_value) max_value = original_point_values[j][selected_vars[i]];
			if (original_point_values[j][selected_vars[i]] < min_value) min_value = original_point_values[j][selected_vars[i]];
		}

		assert(max_value - min_value != 0);

		for (int j = 0; j < original_point_values.size(); ++j) {
			float scale_value = (original_point_values[j][selected_vars[i]] - min_value) / (max_value - min_value);
			X0->set(j, i, scale_value);
		}
	}

	smat::Matrix<double> *D = new smat::Matrix<double>(original_point_values.size(), original_point_values.size(), 0.0);
	float total_weight = 0;
	for (int i = 0; i < selected_vars.size(); ++i) total_weight += var_weights[selected_vars[i]];
	for (int i = 0; i < original_point_values.size() - 1; ++i){
		D->set(i, i, 0);
		for (int j = i + 1; j < original_point_values.size(); ++j) {
			float dis = 0;
			for (int k = 0; k < selected_vars.size(); ++k)
				dis += pow((X0->get(i, k) - X0->get(j, k)) * var_weights[selected_vars[k]], 2) + 1e-20;
			dis = sqrt(dis / total_weight);
			assert(dis != 0);
			D->set(i, j, dis);
			D->set(j, i, dis);
		}
	}

	int dim = 2;
	
	//int iteration = 40;
	int iteration = 40;

	smat::Matrix<double> * X1 = MDS_SMACOF(D, NULL, dim, iteration); // without initialization

	original_point_pos.resize(original_point_values.size());
	for (int i = 0; i < original_point_values.size(); ++i) {
		original_point_pos[i].resize(2);
		original_point_pos[i][0] = X1->get(i, 0) * 10000;
		original_point_pos[i][1] = X1->get(i, 1) * 10000;
	}
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

	max_pos_range = -1e10;
	for (int i = 0; i < vec[0].size(); ++i){
		float minValue = 1e10;
		float maxValue = -1e10;

		for (int j = 0; j < vec.size(); ++j){
			if (minValue > vec[j][i]) minValue = vec[j][i];
			if (maxValue < vec[j][i]) maxValue = vec[j][i];
		}
		if (maxValue - minValue > max_pos_range) max_pos_range = maxValue - minValue;

		ranges[i].resize(2);
		ranges[i][0] = minValue;
		ranges[i][1] = maxValue;
	}

	for (int i = 0; i < vec[0].size(); ++i){
		float mid_value = (ranges[i][0] + ranges[i][1]) / 2;
		
		ranges[i][0] = mid_value - max_pos_range / 2;
		ranges[i][1] = mid_value + max_pos_range / 2;

		for (int j = 0; j < vec.size(); ++j)
			vec[j][i] = (vec[j][i] - mid_value) / max_pos_range + 0.5;
	}
}