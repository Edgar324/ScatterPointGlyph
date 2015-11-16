#include "similarity_extractor.h"
#include <queue>

SimilarityExtractor::SimilarityExtractor() {

}

SimilarityExtractor::~SimilarityExtractor() {

}

void SimilarityExtractor::ExtractCosts(float thres) {
	ExtractProposalGestalt(thres);

	// extract cost
	int label_num = this->proposal_gestalt.size() + 1;
	int gestalt_num = this->proposal_gestalt.size();
	int site_num = this->proposal_gestalt.size();

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
	for (int i = 0; i < gestalt_num; ++i) {
		float average_value = 0;
		int point_count_sum = 0;
		for (int j = 0; j < proposal_gestalt[i].size(); ++j) {
			int site_index = proposal_gestalt[i][j];
			average_value += gestalt_candidates->site_average_value[site_index] * gestalt_candidates->site_point_num[site_index];
			point_count_sum += gestalt_candidates->site_point_num[site_index];
		}
		average_value /= point_count_sum;

		float value_var = 0;
		for (int j = 0; j < proposal_gestalt[i].size(); ++j) {
			int site_index = proposal_gestalt[i][j];
			value_var += pow(gestalt_candidates->site_average_value[i] - average_value, 2) * gestalt_candidates->site_point_num[site_index];

			float value_dis = abs(gestalt_candidates->site_average_value[site_index] - average_value);
			data_cost[site_index][i] = value_dis;
		}
		value_var = sqrt(value_var / point_count_sum);
		label_cost[i] = value_var;
	}
	label_cost[label_num - 1] = 0;
	NormalizeVec(label_cost);

	for (int i = 0; i < site_num; ++i)
		if (gestalt_candidates->is_site_labeled[i]) data_cost[i][label_num - 1] = 0;
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

void SimilarityExtractor::ExtractProposalGestalt(float thres) {
	std::vector< bool > is_selected;

	int gestalt_num = gestalt_candidates->gestalt_candidates.size();

	for (int i = 0; i < gestalt_num; ++i) {
		if (gestalt_candidates->is_site_labeled[i]) continue;

		is_selected.resize(gestalt_candidates->gestalt_candidates[i].size());
		is_selected.assign(is_selected.size(), false);

		int center_index = gestalt_candidates->gestalt_candidates[i][0];

		std::queue< int > extending_queue;
		extending_queue.push(0);
		while (!extending_queue.empty()) {
			int site_index = gestalt_candidates->gestalt_candidates[i][extending_queue.front()];
			is_selected[extending_queue.front()] = true;
			extending_queue.pop();

			for (int j = 0; j < gestalt_candidates->gestalt_candidates[i].size(); ++j)
				if (!is_selected[j]) {
					int next_site_index = gestalt_candidates->gestalt_candidates[i][j];
					if (gestalt_candidates->is_site_labeled[next_site_index]) continue;

					if (gestalt_candidates->site_connecting_status[site_index][next_site_index]) {
						float dis = abs(gestalt_candidates->site_average_value[site_index] - gestalt_candidates->site_average_value[center_index]);
						if (dis < thres) extending_queue.push(j);
					}
				}
		}
		this->proposal_gestalt[i].clear();
		for (int j = 0; j < is_selected.size(); ++j)
			if (is_selected[j]) this->proposal_gestalt[i].push_back(gestalt_candidates->gestalt_candidates[i][j]);
	}
}