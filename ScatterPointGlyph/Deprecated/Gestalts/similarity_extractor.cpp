#include "similarity_extractor.h"
#include <queue>

SimilarityExtractor::SimilarityExtractor() : PropertyExtractor() {

}

SimilarityExtractor::~SimilarityExtractor() {

}

void SimilarityExtractor::ExtractCosts(float thres) {
	ExtractProposalGestalt(thres);

	// extract cost
	int label_num = this->proposal_gestalt.size();
	int site_num = gestalt_candidates->site_nodes.size();
	int var_num = gestalt_candidates->var_weights.size();

	label_cost.resize(label_num);
	label_cost.assign(label_num, -1);

	smooth_cost.resize(label_num);
	for (int i = 0; i < label_num; ++i) {
		smooth_cost[i].resize(label_num);
		smooth_cost[i].assign(label_num, -1);
	}

	data_cost.resize(site_num);
	for (int i = 0; i < site_num; ++i) {
		data_cost[i].resize(label_num);
		data_cost[i].assign(label_num, -1);
	}

	for (int i = 0; i < label_num; ++i) {
		std::vector< float > average_value;
		average_value.resize(var_num, 0);

		int point_count_sum = 0;
		for (int j = 0; j < proposal_gestalt[i].size(); ++j) {
			int site_index = proposal_gestalt[i][j];
			for (int k = 0; k < var_num; ++k)
				average_value[k] += gestalt_candidates->site_nodes[site_index]->average_values[k];
			point_count_sum += 1;
		}
		for (int j = 0; j < var_num; ++j) average_value[j] /= point_count_sum;

		float value_var = 0;
		for (int j = 0; j < proposal_gestalt[i].size(); ++j) {
			int site_index = proposal_gestalt[i][j];
			float value_dis = 0;
			for (int k = 0; k < var_num; ++k)
				value_dis += abs(gestalt_candidates->site_nodes[site_index]->average_values[k] - average_value[k]) * gestalt_candidates->var_weights[k];
			value_var += pow(value_dis, 2);

			data_cost[site_index][i] = value_dis;
		}
		value_var = sqrt(value_var / point_count_sum);
		label_cost[i] = value_var;
	}

	for (int i = 0; i < label_num; ++i){
		smooth_cost[i][i] = 0;

		for (int j = i + 1; j < label_num; ++j) {
			smooth_cost[i][j] = 0;
			for (int k = 0; k < var_num; ++k)
				smooth_cost[i][j] += abs(gestalt_candidates->clusters[i]->average_values[k] - gestalt_candidates->clusters[j]->average_values[k]) * gestalt_candidates->var_weights[k];
			smooth_cost[j][i] = smooth_cost[i][j];
		}
	}

	for (int i = 0; i < label_cost.size(); ++i) {
		if (label_cost[i] < 0) label_cost[i] = this->maximum_cost;
		label_cost[i] = label_cost[i] * scales[2] + bias[2];
	}
	for (int i = 0; i < data_cost.size(); ++i)
		for (int j = 0; j < data_cost[i].size(); ++j) {
			if (data_cost[i][j] < 0) data_cost[i][j] = this->maximum_cost;
			data_cost[i][j] = data_cost[i][j] * scales[0] + bias[0];
		}
	for (int i = 0; i < label_num; ++i){
		for (int j = i + 1; j < label_num; ++j) {
			smooth_cost[i][j] = smooth_cost[i][j] * scales[1] + bias[1];
			smooth_cost[j][i] = smooth_cost[i][j];
		}
	}
}

void SimilarityExtractor::ExtractProposalGestalt(float thres) {
	std::vector< bool > is_selected;

	int gestalt_num = gestalt_candidates->gestalt_cluster_index.size();
	this->proposal_gestalt.clear();
	this->proposal_clusters.clear();

	for (int i = 0; i < gestalt_num; ++i) {
		is_selected.resize(gestalt_candidates->gestalt_cluster_index[i].size());
		is_selected.assign(is_selected.size(), false);

		int center_index = gestalt_candidates->gestalt_cluster_index[i][0];

		std::queue< int > extending_queue;
		extending_queue.push(0);
		while (!extending_queue.empty()) {
			int site_index = gestalt_candidates->gestalt_cluster_index[i][extending_queue.front()];
			is_selected[extending_queue.front()] = true;
			extending_queue.pop();

			for (int j = 0; j < gestalt_candidates->gestalt_cluster_index[i].size(); ++j)
				if (!is_selected[j]) {
					int next_site_index = gestalt_candidates->gestalt_cluster_index[i][j];

					if (gestalt_candidates->cluster_connecting_status[site_index][next_site_index]) {
						float dis = 0;
						for (int k = 0; k < gestalt_candidates->var_weights.size(); ++k) 
							dis += abs(gestalt_candidates->clusters[site_index]->average_values[k] - gestalt_candidates->gestalt_average_values[center_index][k]) *  gestalt_candidates->var_weights[k];
						if (dis < thres) extending_queue.push(j);
					}
				}
		}

		std::vector< int > proposal;
		std::vector< int > cluster_proposal;
		for (int j = 0; j < is_selected.size(); ++j)
			if (is_selected[j]) {
				cluster_proposal.push_back(gestalt_candidates->gestalt_cluster_index[i][j]);
				for (int k = 0; k < gestalt_candidates->site_nodes.size(); ++k)
					if (gestalt_candidates->basic_node_index[k] == gestalt_candidates->gestalt_cluster_index[i][j]) proposal.push_back(k);
			}
		this->proposal_clusters.push_back(cluster_proposal);
		this->proposal_gestalt.push_back(proposal);
	}
}