/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef PARALLEL_DATASET_H_
#define PARALLEL_DATASET_H_

#include "view_dataset.h"
#include <vector>
using namespace std;
#include <QtGui/QColor>

class ParallelRecord 
{
public:
    ParallelRecord() { }        
    ~ParallelRecord() { }

    vector<float> values;
};

class ParallelDataset : public ViewDataset
{
    Q_OBJECT

public:
    ParallelDataset();
    ~ParallelDataset();

    bool CompleteInput();
	void UpdateGaussian();
    virtual void Clear();
    virtual void Modified();

    // attributes which must be set
    vector<QString > subset_names;
    vector<vector<ParallelRecord*>> subset_records;
    vector<vector<QString>> axis_anchors;
    vector<QString > axis_names;
	vector<QColor > subset_colors;


    // attributes which can be set automatically
	vector<vector<float>> var_centers;
	vector<vector<float>> var_width;
	vector<vector<float>> var_std_dev;
    vector<vector<QColor>> record_color;
    vector<vector<bool>> is_record_selected;

    vector<bool> is_subset_visible;
    vector<float> subset_opacity;
    vector<bool> is_axis_selected;
    vector<int> mapped_axis;

	bool is_axis_weight_enabled;
	vector<float> axis_weights;

    bool is_cluster_enabled;
    vector<vector<ParallelRecord*>> cluster_centers;

    bool is_range_filter_enabled;
    vector<vector<float>> axis_value_filter_range;

    bool is_edge_bundling_enabled;
    bool is_correlation_analysis_enabled;
	bool is_gaussian_enabled;

	bool is_updating;
};

#endif