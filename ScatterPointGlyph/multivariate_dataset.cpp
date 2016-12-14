/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "multivariate_dataset.h"
#include <fstream>
#include <QStringList>
#include "data_projector.h"
#include "color_mapping_generator.h"

MultivariateDataset::MultivariateDataset(const char* file_name) {
    is_good_ = LoadFromFile(file_name);
}

MultivariateDataset::~MultivariateDataset() {

}

double* MultivariateDataset::GetRecordPos(int i) {
    temp_pos_[0] = projected_pos_(i, 0);
    temp_pos_[1] = projected_pos_(i, 1);
    return temp_pos_;
}

void MultivariateDataset::GetPosRanges(double& left, double& right, double& bottom, double& top) {
    left = 0;
    right = 1;
    bottom = 0;
    top = 1;
}

void MultivariateDataset::GetAveragePos(const vector<int>& record_ids, double& x, double& y) {
    x = 0;
    y = 0;
	for (int i = 0; i < record_ids.size(); ++i) {
        x += projected_pos_(record_ids[i], 0);
        y += projected_pos_(record_ids[i], 1);
	}
    if (record_ids.size() != 0) {
        x /= record_ids.size();
        y /= record_ids.size();
    }
}

void MultivariateDataset::GetAverageValues(const vector<int>& record_ids, vector<double>& values, vector<double>& std_deviations) {
    values.resize(var_num_);
    values.assign(var_num_, 0);
    std_deviations.resize(var_num_);
    std_deviations.assign(var_num_, 0);

    for (int i = 0; i < record_ids.size(); ++i)
		for (int j = 0; j < var_num_; ++j)
	        values[j] += records_(j, record_ids[i]);
	for (int i = 0; i < var_num_; ++i)
		values[i] /= record_ids.size();

	for (int i = 0; i < record_ids.size(); ++i)
		for (int j = 0; j < var_num_; ++j)
			std_deviations[j] += pow(records_(j, record_ids[i]) - values[j], 2);

	for (int i = 0; i < var_num_; ++i)
        if (record_ids.size() > 1)
            std_deviations[i] = sqrt(std_deviations[i] / (record_ids.size() - 1));
        else
		    std_deviations[i] = sqrt(std_deviations[i] / record_ids.size());
}

double MultivariateDataset::Distance(int record_one, int record_two, double ratio /*= 1.0*/) {
    return 0.0;
}

bool MultivariateDataset::ApplyProjection(DataProjector* projector) {
    projector->Project(records_, projected_pos_);
    is_projected_ = true;
    
    // Normalize the projections
    pos_ranges_.resize(2, 2);
    pos_ranges_(0, 0) = 1e10;
    pos_ranges_(0, 1) = -1e10;
    pos_ranges_(1, 0) = 1e10;
    pos_ranges_(1, 1) = -1e10;

    for (int i = 0; i < record_num_; i++) {
        if (projected_pos_(i, 0) < pos_ranges_(0, 0)) pos_ranges_(0, 0) = projected_pos_(i, 0);
        if (projected_pos_(i, 0) > pos_ranges_(0, 1)) pos_ranges_(0, 1) = projected_pos_(i, 0);

        if (projected_pos_(i, 1) < pos_ranges_(1, 0)) pos_ranges_(1, 0) = projected_pos_(i, 1);
        if (projected_pos_(i, 1) > pos_ranges_(1, 1)) pos_ranges_(1, 1) = projected_pos_(i, 1);
    }

    double max_range = pos_ranges_(0, 1) - pos_ranges_(0, 0);
    if (pos_ranges_(1, 1) - pos_ranges_(1, 0) > max_range)
        max_range = pos_ranges_(1, 1) - pos_ranges_(1, 0);
    double center_x = (pos_ranges_(0, 1) + pos_ranges_(0, 0)) / 2;
    double center_y = (pos_ranges_(1, 1) + pos_ranges_(1, 0)) / 2;

    for (int i = 0; i < record_num_; i++) {
        projected_pos_(i, 0) = (projected_pos_(i, 0) - center_x) / max_range;
        projected_pos_(i, 1) = (projected_pos_(i, 1) - center_y) / max_range;
    }

    pos_ranges_(0, 0) = center_x - max_range / 2;
    pos_ranges_(0, 1) = center_x + max_range / 2;
    pos_ranges_(1, 0) = center_y - max_range / 2;
    pos_ranges_(1, 1) = center_y + max_range / 2;

    return true;
}

bool MultivariateDataset::LoadFromFile(const char* file_name) {
    ifstream input_file(file_name);

    if (!input_file.good()) return false;

	char char_str[1000];
	input_file.getline(char_str, 1000);
	QString value_str = QString::fromLocal8Bit(char_str);
	QStringList value_list = value_str.split(' ');

	record_num_ = value_list.at(0).toInt();
	var_num_ = value_list.at(1).toInt();

	input_file.getline(char_str, 1000);
	value_str = QString::fromLocal8Bit(char_str);
	value_list = value_str.split(' ');
	for (int i = 0; i < value_list.size(); ++i) 
        var_names_.push_back(value_list.at(i));

    var_weights_.resize(var_num_);
    var_weights_.assign(var_num_, 1.0f / var_num_);

    ColorMappingGenerator::GetInstance()->GetQualitativeColors(var_num_, var_colors_);

    records_.resize(var_num_, record_num_);

	for (int i = 0; i < record_num_; ++i) {
		input_file.getline(char_str, 1000);
		value_str = QString::fromLocal8Bit(char_str);
		value_list = value_str.split(',');

		for (int j = 0; j < var_num_; ++j)
			records_(j, i) = value_list.at(j).toFloat();
	}

    var_ranges_.resize(var_num_, 2);

	for (int i = 0; i < var_num_; ++i){
		float min_value = 1e10;
		float max_value = -1e10;

		for (int j = 0; j < record_num_; ++j){
			if (min_value > records_(i, j)) min_value = records_(i, j);
			if (max_value < records_(i, j)) max_value = records_(i, j);
		}

		var_ranges_(i, 0) = min_value;
		var_ranges_(i, 1) = max_value;
	}

    // normalize the data
    for (int i = 0; i < var_num_; i++) {
        double temp_range = var_ranges_(i, 1) - var_ranges_(i, 0) + 1e-10;
        for (int j = 0; j < record_num_; j++) {
            records_(i, j) = (records_(i, j) - var_ranges_(i, 0)) / temp_range;
        }
    }

	input_file.close();

    return true;
}