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
	for (int i = 0; i < label_num; ++i) {
		float average_edge_length = 0;
		int edge_count = 0;
		for (int j = 0; j < proposal_gestalt[i].size() - 1; ++j)
			for (int k = j + 1; k < proposal_gestalt[i].size(); ++k) {
				int site_index = proposal_gestalt[i][j];
				int next_node_index = proposal_gestalt[i][k];
				if (gestalt_candidates->site_connecting_status[site_index][next_node_index]) {
					average_edge_length += sqrt(pow(gestalt_candidates->site_nodes[site_index]->center_pos[0] - gestalt_candidates->site_nodes[next_node_index]->center_pos[0], 2)
						+ pow(gestalt_candidates->site_nodes[site_index]->center_pos[1] - gestalt_candidates->site_nodes[next_node_index]->center_pos[1], 2));
					edge_count++;
				}
			}
		if (edge_count != 0) average_edge_length /= edge_count;
		label_cost[i] = average_edge_length + 1e-3;
	}
	NormalizeVec(label_cost);

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
	NormalizeVec(data_cost);

	// update smooth cost
	for (int i = 0; i < label_num; ++i){
		smooth_cost[i][i] = 0;

		for (int j = i + 1; j < label_num; ++j) {
			smooth_cost[i][j] = 1.0;
			smooth_cost[j][i] = smooth_cost[i][j];
		}
	}
}

void ProximityExtractor::ExtractProposalGestalt(float thres) {
	std::vector< bool > is_selected;

	int gestalt_num = gestalt_candidates->gestalt_cluster_index.size();
	this->proposal_gestalt.clear();

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
						float dis = sqrt(pow(gestalt_candidates->clusters[site_index]->center_pos[0] - gestalt_candidates->clusters[next_site_index]->center_pos[0], 2)
							+ pow(gestalt_candidates->clusters[site_index]->center_pos[1] - gestalt_candidates->clusters[next_site_index]->center_pos[1], 2));
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