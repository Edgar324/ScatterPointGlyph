/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "parallel_dataset_builder.h"
#include "color_mapping_generator.h"
#include "parallel_dataset.h"
#include "scatter_point_dataset.h"

void ParallelDatasetBuilder::Build(ScatterPointDataset* point_dataset, vector<int>& selected_var_index, 
    vector<vector<int>>& point_ids, ParallelDataset* parallel_dataset) {
    parallel_dataset->Clear();

    int var_num = selected_var_index.size();
    int cluster_num = point_ids.size();
	parallel_dataset->subset_names.resize(cluster_num);
	parallel_dataset->subset_colors.resize(cluster_num);
	parallel_dataset->subset_records.resize(cluster_num);

	parallel_dataset->axis_names.resize(var_num);
	parallel_dataset->axis_anchors.resize(var_num);

    vector<QColor> cluster_colors;
    ColorMappingGenerator::GetInstance()->GetQualitativeColors(cluster_num, cluster_colors);

	for (int i = 0; i < point_ids.size(); ++i) {
		parallel_dataset->subset_names[i] = QString("Cluster %0").arg(i);
        parallel_dataset->subset_colors[i] = cluster_colors[i];
        for (int j = 0; j < point_ids[i].size(); ++j) {
            ParallelRecord* record = new ParallelRecord;
            record->values.resize(var_num);
            for (int k = 0; k < var_num; ++k)
                record->values[k] = point_dataset->normalized_point_values[selected_var_index[k]][point_ids[i][j]];
			parallel_dataset->subset_records[i].push_back(record);
        }
	}

	for (int i = 0; i < var_num; ++i) {
		parallel_dataset->axis_names[i] = point_dataset->var_names[selected_var_index[i]];
		parallel_dataset->axis_anchors[i].push_back(QString("%0").arg(point_dataset->original_value_ranges[selected_var_index[i]][0]));
		parallel_dataset->axis_anchors[i].push_back(QString("%0").arg(point_dataset->original_value_ranges[selected_var_index[i]][1]));
	}
	parallel_dataset->CompleteInput();
	parallel_dataset->UpdateGaussian();
}