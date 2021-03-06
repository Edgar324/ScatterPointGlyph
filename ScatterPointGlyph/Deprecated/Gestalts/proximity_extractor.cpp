#include "proximity_extractor.h"
#include <queue>

ProximityExtractor::ProximityExtractor()
	: PropertyExtractor() {

}

ProximityExtractor::~ProximityExtractor() {

}

void ProximityExtractor::ExtractCosts(float thres) {
	// extract proposal gestalt
	ExtractProposalGestalt(thres);

	// extract cost
	int label_num = this->proposal_gestalt.size();
	int site_num = gestalt_candidates->site_nodes.size();

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

	// update label cost
	std::vector< float > average_edge_length;
	average_edge_length.resize(label_num, 0);
	for (int i = 0; i < label_num; ++i) {
		average_edge_length[i] = 0;
		int edge_count = 0;
		for (int j = 0; j < proposal_gestalt[i].size() - 1; ++j)
			for (int k = j + 1; k < proposal_gestalt[i].size(); ++k) {
				int site_index = proposal_gestalt[i][j];
				int next_node_index = proposal_gestalt[i][k];
				if (gestalt_candidates->site_connecting_status[site_index][next_node_index]) {
					average_edge_length[i] += sqrt(pow(gestalt_candidates->site_nodes[site_index]->center_pos[0] - gestalt_candidates->site_nodes[next_node_index]->center_pos[0], 2)
						+ pow(gestalt_candidates->site_nodes[site_index]->center_pos[1] - gestalt_candidates->site_nodes[next_node_index]->center_pos[1], 2));
					edge_count++;
				}
			}
		if (edge_count != 0) average_edge_length[i] /= edge_count;
		label_cost[i] = average_edge_length[i];
	}

	// update data cost
	std::vector< std::vector< float > > label_center;
	label_center.resize(label_num);
	for (int i = 0; i < label_num; ++i) {
		label_center[i].resize(2);
		label_center[i][0] = 0;
		label_center[i][1] = 0;
		int point_count_sum = 0;

		for (int j = 0; j < proposal_gestalt[i].size(); ++j) {
			int site_index = proposal_gestalt[i][j];
			label_center[i][0] += gestalt_candidates->site_nodes[site_index]->center_pos[0];
			label_center[i][1] += gestalt_candidates->site_nodes[site_index]->center_pos[1];
			point_count_sum += 1;
		}

		label_center[i][0] /= point_count_sum;
		label_center[i][1] /= point_count_sum;
	}

	for (int i = 0; i < label_num; ++i) {
		for (int j = 0; j < proposal_gestalt[i].size(); ++j) {
			int site_index = proposal_gestalt[i][j];
			float dis = sqrt(pow(gestalt_candidates->site_nodes[site_index]->center_pos[0] - label_center[i][0], 2)
				+ pow(gestalt_candidates->site_nodes[site_index]->center_pos[1] - label_center[i][1], 2));
			data_cost[site_index][i] = dis + 1e-3;
		}
	}

	for (int i = 0; i < label_num; ++i){
		smooth_cost[i][i] = 0;

		for (int j = i + 1; j < label_num; ++j) {
			smooth_cost[i][j] = 0.0;
			smooth_cost[i][j] = abs(average_edge_length[i] - average_edge_length[j]);
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

void ProximityExtractor::ExtractProposalGestalt(float thres) {
	std::vector< bool > is_selected;

	int gestalt_num = gestalt_candidates->gestalt_cluster_index.size();
	this->proposal_gestalt.clear();
	this->proposal_clusters.clear();

	std::vector< std::vector< float > > min_dis_vec;
	min_dis_vec.resize(gestalt_candidates->clusters.size());
	for (int i = 0; i < gestalt_candidates->clusters.size(); ++i) min_dis_vec[i].resize(gestalt_candidates->clusters.size(), 1e10);
	for (int i = 0; i < gestalt_candidates->basic_node_index.size() - 1; ++i)
		for (int j = i + 1; j < gestalt_candidates->basic_node_index.size(); ++j) {
			if (!gestalt_candidates->site_connecting_status[i][j]) continue;
			int site_index = gestalt_candidates->basic_node_index[i];
			int next_site_index = gestalt_candidates->basic_node_index[j];
			float temp_dis = sqrt(pow(gestalt_candidates->site_nodes[i]->center_pos[0] - gestalt_candidates->site_nodes[j]->center_pos[0], 2)
				+ pow(gestalt_candidates->site_nodes[i]->center_pos[1] - gestalt_candidates->site_nodes[j]->center_pos[1], 2));
			if (min_dis_vec[site_index][next_site_index] > temp_dis) {
				min_dis_vec[site_index][next_site_index] = temp_dis;
				min_dis_vec[next_site_index][site_index] = temp_dis;
			}
		}

	for (int i = 0; i < gestalt_num; ++i) {
		is_selected.resize(gestalt_candidates->gestalt_cluster_index[i].size());
		is_selected.assign(is_selected.size(), false);
		
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
						/*float dis = sqrt(pow(gestalt_candidates->clusters[site_index]->center_pos[0] - gestalt_candidates->clusters[next_site_index]->center_pos[0], 2)
							+ pow(gestalt_candidates->clusters[site_index]->center_pos[1] - gestalt_candidates->clusters[next_site_index]->center_pos[1], 2));*/
						if (min_dis_vec[site_index][next_site_index] < thres) extending_queue.push(j);
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