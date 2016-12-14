#include "multi_label_processor.h"
#include <time.h>
#include <queue>
#include <map>
#include <iostream>
#include "./gco-v3.0/GCoptimization.h"
#include "./gco-v3.0/energy.h"

MultiLabelProcessor::MultiLabelProcessor(double label_cost_rate/* = 1.0*/, double data_dis_scale/* = 0.5*/)
    : label_cost_rate_(label_cost_rate), data_dis_scale_(data_dis_scale) {

}

MultiLabelProcessor::~MultiLabelProcessor() {

}

void MultiLabelProcessor::ExtractEstimatedModels(const vector<vector<bool>>& edges) {
    this->estimated_models_.clear();

    int point_num = dis_mat_.size();

	vector<bool> is_reached;
	is_reached.resize(point_num);

	for (int i = 0; i < point_num; ++i) {
		is_reached.assign(point_num, false);

		dis_mat_[i][i] = 0;
        is_reached[i] = true;

		for (int j = 0; j < point_num; ++j) {
			float min_dist = 1e10;
			int min_index = -1;

			for (int k = 0; k < point_num; ++k)
				if (!is_reached[k] && dis_mat_[i][k] < min_dist && dis_mat_[i][k] < max_radius_) {
					min_dist = dis_mat_[i][k];
					min_index = k;
				}
			if (min_index == -1) break;
			is_reached[min_index] = true;

			for (int k = 0; k < point_num; ++k)
				if (!is_reached[k] && edges[min_index][k]) {
					if (dis_mat_[i][k] > dis_mat_[i][min_index] + dis_mat_[min_index][k]) dis_mat_[i][k] = dis_mat_[i][min_index] + dis_mat_[min_index][k];
				}
		}

		vector<int> model;
		for (int j = 0; j < is_reached.size(); ++j)
			if (is_reached[j]) model.push_back(j);
		this->estimated_models_.push_back(model);
	}
}

void MultiLabelProcessor::GenerateCluster(const vector<vector<bool>>& edges, vector<vector<int>>& clusters) {
	try {
		int point_num = edges.size();
		int site_num = edges.size();

		// Step 1: construct class
		GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(site_num, point_num);

		// Step 2: multi-label graph cut
		for (int i = 0; i < site_num - 1; ++i)
			for (int j = i + 1; j < site_num; ++j) {
				if (edges[i][j]) gc->setNeighbors(i, j);
			}

		for (int i = 0; i < site_num; ++i) {
			for (int j = 0; j < point_num; ++j) {
				gc->setDataCost(i, j, (int)(data_cost_[i][j] * 1000));
			}
		}
		for (int i = 0; i < point_num; ++i)
			for (int j = 0; j < point_num; ++j) {
				gc->setSmoothCost(i, j, (int)(smooth_cost_[i][j] * 1000));
			}

		vector<int> temp_label_cost;
		temp_label_cost.resize(point_num);
		for (int i = 0; i < point_num; ++i) temp_label_cost[i] = (int)(label_cost_[i] * 1000);
		gc->setLabelCost(temp_label_cost.data());

		gc->expansion(2); // run expansion for 2 iterations. For swap use gc->swap(num_iterations);

        cout << "Optimizatoin finished!" << endl;

        vector<int> result_label(site_num);
		for (int i = 0; i < site_num; ++i) result_label[i] = gc->whatLabel(i);

        // Step 3: adapt results for the clusters
        map<int, int> label_map;
        for (int i = 0; i < result_label.size(); i++) {
            if (label_map.find(result_label[i]) == label_map.end()) {
                label_map.insert(map<int, int>::value_type(result_label[i], (int)label_map.size()));
            }
        }
        clusters.clear();
        clusters.resize(label_map.size());
        for (int i = 0; i < result_label.size(); i++)
            clusters[label_map[result_label[i]]].push_back(i);

		delete gc;
	}
	catch (GCException e) {
		e.Report();
	}
}

void MultiLabelProcessor::GenerateCluster(const vector<vector<double>>& pos, const vector<vector<double>>& value, 
    const vector<double>& weights, const vector<vector<bool>>& edges, 
    double radius, vector<vector<int>>& clusters) {

    int point_num = pos.size();
    max_radius_ = radius;

    // Update distance matrix for model extraction
    dis_mat_.clear();
    dis_mat_.resize(point_num, vector<double>(point_num, 1e10));
    average_value_dis_ = 0.0;
    int connected_edge_num = 0;
    double pos_dis, value_dis;
	for (int i = 0; i < point_num - 1; ++i)
		for (int j = i + 1; j < point_num; ++j)
            if (edges[i][j]) {
                pos_dis = sqrt(pow(pos[i][0] - pos[j][0], 2) + pow(pos[i][1] - pos[j][1], 2));
                value_dis = 0;
                for (int k = 0; k < weights.size(); ++k)
                    value_dis += abs(value[i][k] - value[j][k]) * weights[k];
                average_value_dis_ += value_dis;
                connected_edge_num++;

                dis_mat_[i][j] = data_dis_scale_ * value_dis + (1.0 - data_dis_scale_) * pos_dis;
                dis_mat_[j][i] = dis_mat_[i][j];
            }
    if (connected_edge_num != 0) average_value_dis_ /= connected_edge_num;

    // Extract models using the distance matrix
    ExtractEstimatedModels(edges);

    // Generate energies
	label_cost_.resize(point_num);
	label_cost_.assign(point_num, -1);

	smooth_cost_.resize(point_num);
	for (int i = 0; i < point_num; ++i) {
		smooth_cost_[i].resize(point_num);
		smooth_cost_[i].assign(point_num, -1);
	}

	data_cost_.resize(point_num);
	for (int i = 0; i < point_num; ++i) {
		data_cost_[i].resize(point_num);
		data_cost_[i].assign(point_num, -1);
	}

	vector<vector<double>> average_values;
	average_values.resize(point_num);

	for (int i = 0; i < point_num; ++i) {
		int model_size = estimated_models_[i].size();

		average_values[i].resize(weights.size(), 0);

		float center_x = 0, center_y = 0;
		for (int j = 0; j < model_size; ++j) {
			int site_index = estimated_models_[i][j];
			for (int k = 0; k < weights.size(); ++k)
				average_values[i][k] += value[site_index][k];
			center_x += pos[site_index][0];
			center_y += pos[site_index][1];
		}
		for (int k = 0; k < weights.size(); ++k)
			average_values[i][k] /= model_size;
		center_x /= model_size;
		center_y /= model_size;

		float value_var = 0;
		for (int j = 0; j < model_size; ++j) {
			int site_index = estimated_models_[i][j];
			double value_dis = 0;
			for (int k = 0; k < weights.size(); ++k)
				value_dis += abs(value[site_index][k] - average_values[i][k]) * weights[k];
			double pos_dis = sqrt(pow(pos[site_index][0] - center_x, 2) + pow(pos[site_index][1] - center_y, 2));

			data_cost_[site_index][i] = data_dis_scale_ * value_dis + (1.0 - data_dis_scale_) * pos_dis;

			value_var += pow(value_dis, 2);
		}
        if (model_size > 1)
            value_var = sqrt(value_var / (model_size - 1));
        else
		    value_var = sqrt(value_var / model_size);
		label_cost_[i] = label_cost_rate_* value_var * point_num + 0.01;
	}

	for (int i = 0; i < point_num; ++i){
		smooth_cost_[i][i] = 0;

		for (int j = i + 1; j < point_num; ++j) {
            smooth_cost_[i][j] = average_value_dis_;

			smooth_cost_[j][i] = smooth_cost_[i][j];
		}
	}

	for (int i = 0; i < label_cost_.size(); ++i) {
		if (label_cost_[i] < 0) label_cost_[i] = 100.0;
	}
	for (int i = 0; i < data_cost_.size(); ++i)
		for (int j = 0; j < data_cost_[i].size(); ++j) {
			if (data_cost_[i][j] < 0) data_cost_[i][j] = 100.0;
		}

    GenerateCluster(edges, clusters);
}

void MultiLabelProcessor::GenerateCluster(const vector<vector<double>>& value, 
    const vector<double>& weights, const vector<vector<bool>>& edges, 
    double radius, vector<vector<int>>& clusters) {
    cout << "Currently not implemented!" << endl;
}
