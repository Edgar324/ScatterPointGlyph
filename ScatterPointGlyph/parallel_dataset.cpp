/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "parallel_dataset.h"

ParallelDataset::ParallelDataset()
	: is_edge_bundling_enabled(false), is_correlation_analysis_enabled(false),
	is_cluster_enabled(false), is_range_filter_enabled(false), is_axis_weight_enabled(false),
	is_gaussian_enabled(false), is_updating(false) {
}

ParallelDataset::~ParallelDataset(){

}

bool ParallelDataset::CompleteInput(){
    if ( subset_names.size() != subset_records.size() || axis_anchors.size() != axis_names.size() ) return false;

    if ( subset_colors.size() != subset_names.size() ){
        subset_colors.clear();
        /// TODO:select colors here
		for (int i = 0; i < subset_names.size(); ++i)
			subset_colors.push_back(QColor(200, 200, 200));
    }

	if (subset_colors.size() == subset_names.size()) {
		record_color.resize(subset_names.size());
		for (int i = 0; i < subset_names.size(); ++i) {
			record_color[i].resize(subset_records[i].size());
			record_color[i].assign(subset_records[i].size(), subset_colors[i]);
		}
	}

    if ( is_subset_visible.size() != subset_names.size() ){
        is_subset_visible.clear();
        is_subset_visible.resize(subset_names.size(), true);
    }

    if ( subset_opacity.size() != subset_names.size() ){
        subset_opacity.clear();
        subset_opacity.resize(subset_names.size(), 1.0);
    }

    if ( mapped_axis.size() != axis_names.size() ){
        mapped_axis.resize(axis_names.size());
        for ( int i = 0; i < mapped_axis.size(); ++i ) mapped_axis[i] = i;
    }

    if ( is_record_selected.size() != subset_names.size() ){
        is_record_selected.resize(subset_names.size());
        for ( int i = 0; i < is_record_selected.size(); ++i )
            if ( is_record_selected[i].size() != subset_records[i].size() ){
                is_record_selected[i].resize(subset_records[i].size());
                is_record_selected[i].assign(is_record_selected[i].size(), true);
                //memset(&is_record_selected[i][0], 0, subset_records[i].size());
            }
    }

    if ( is_axis_selected.size() != axis_names.size() ){
        is_axis_selected.clear();
        is_axis_selected.resize(axis_names.size(), false);
    }
    return true;
}

void ParallelDataset::UpdateGaussian() {
	// gausian analysis
	var_centers.resize(subset_names.size());
	var_width.resize(subset_names.size());
	var_std_dev.resize(subset_names.size());
	for (int i = 0; i < subset_names.size(); ++i) {
		var_centers[i].resize(axis_names.size());
		var_centers[i].assign(axis_names.size(), 0);
		var_width[i].resize(axis_names.size());
		var_width[i].assign(axis_names.size(), 0);
		var_std_dev[i].resize(axis_names.size());
		var_std_dev[i].assign(axis_names.size(), 0);

        if (this->subset_records[i].size() == 0) continue;

		for (int j = 0; j < axis_names.size(); ++j) {
			float average = 0;
			std::vector<float> sub_value;
			sub_value.resize(this->subset_records[i].size(), 1000);
			int accu_count = 0;
			for (int k = 0; k < this->subset_records[i].size(); ++k) {
				if (!this->is_record_selected[i][k]) continue;
				sub_value[k] = subset_records[i][k]->values[j];
				average += subset_records[i][k]->values[j];
				accu_count++;
			}
			average /= accu_count;
			for (int k = 0; k < this->subset_records[i].size(); ++k)
				if (this->is_record_selected[i][k]) sub_value[k] = abs(sub_value[k] - average);
			std::sort(sub_value.begin(), sub_value.end());
			var_centers[i][j] = average;
			var_width[i][j] = sub_value[(int)(0.8 * (accu_count - 1))];
			if (var_width[i][j] < 0.01) var_width[i][j] = 0.01;
			for (int k = 0; k < sub_value.size(); ++k)
				var_std_dev[i][j] += pow(sub_value[k], 2);
			var_std_dev[i][j] = sqrt(var_std_dev[i][j] / sub_value.size());
			if (var_std_dev[i][j] < 0.01) var_std_dev[i][j] = 0.01;
		}
	}

	//if (subset_names.size() >= 1) is_gaussian_enabled = true;
}

void ParallelDataset::Clear(){
    this->subset_names.clear();
    for ( int i = 0; i < subset_records.size(); ++i )
        for ( int j = 0; j < subset_records[i].size(); ++j ) delete subset_records[i][j];
    subset_records.clear();
    axis_anchors.clear();
    axis_names.clear();
	subset_colors.clear();
	var_centers.clear();
	var_width.clear();

	record_color.clear();
	is_record_selected.clear();
	subset_opacity.clear();
	is_axis_selected.clear();
}

void ParallelDataset::Modified() {
    emit DataChanged();
}