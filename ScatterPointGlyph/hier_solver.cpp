#include "hier_solver.h"
#include <time.h>

HierSolver::HierSolver() 
	: current_cluster_num_(-1), cluster_num_(-1) {

}

HierSolver::~HierSolver() {

}

void HierSolver::SetData(std::vector< std::vector< float > >& value, std::vector< float >& weight) {
	value_ = value;
	weight_ = weight;

	cluster_index_.resize(value_.size());
	for (size_t i = 0; i < value_.size(); ++i) cluster_index_[i] = i;

	pre_cluster_index_[0] = -1;
	pre_cluster_index_[1] = -1;
	current_cluster_num_ = -1;
}

void HierSolver::SetInitialClusterNum(int num){
	cluster_num_ = num;
	current_cluster_num_ = cluster_num_;

	this->KmeansCluster();
}

void HierSolver::SetClusterNum(int num) {
	cluster_num_ = num;
}

void HierSolver::GetPreClusterIndex(int& cluster_one, int& cluster_two){
	cluster_one = pre_cluster_index_[0];
	cluster_two = pre_cluster_index_[1];
}

void HierSolver::GetGlyphData(std::vector< float >& glyph_pos, std::vector< std::vector< float > >& glyph_values){
	glyph_pos.resize(cluster_center_.size() * 2);
	glyph_values.resize(cluster_center_.size());
	for (int i = 0; i < cluster_center_.size(); ++i) {
		glyph_pos[i * 2] = cluster_center_[i][0];
		glyph_pos[i * 2 + 1] = cluster_center_[i][1];

		glyph_values[i].resize(1);
		glyph_values[i][0] = cluster_node_count_[i];
	}
}

void HierSolver::run() {
	pre_cluster_index_[0] = -1;
	pre_cluster_index_[1] = -1;
	if (current_cluster_num_ <= cluster_num_) return;

	float min_dis = 1e10;
	for (int i = 0; i < cluster_center_.size() - 1; ++i)
		for (int j = i + 1; j < cluster_center_.size(); ++j) {
			float temp_dis = 0;
			for (int k = 0; k < weight_.size(); ++k)
				temp_dis += abs(cluster_center_[i][k] - cluster_center_[j][k]) * weight_[k];
			if (temp_dis < min_dis) {
				min_dis = temp_dis;
				pre_cluster_index_[0] = i;
				pre_cluster_index_[1] = j;
			}
		}
	
	emit CombinedClusterChanged(pre_cluster_index_[0], pre_cluster_index_[1]);
	this->msleep(1000);

	std::vector< float > temp_new_center;
	temp_new_center.resize(weight_.size(), 0);
	int total_node_count = 0;
	for (int i = 0; i < value_.size(); ++i)
		if (cluster_index_[i] == pre_cluster_index_[0] || cluster_index_[i] == pre_cluster_index_[1]) {
			for (int j = 0; j < weight_.size(); ++j) temp_new_center[j] += value_[i][j];
			total_node_count++;
		}
	if (total_node_count != 0) {
		for (int i = 0; i < weight_.size(); ++i) temp_new_center[i] /= total_node_count;
	}

	// update center and cluster index
	cluster_center_[pre_cluster_index_[0]] = temp_new_center;
	for (int i = 0; i < value_.size(); ++i) {
		if (cluster_index_[i] == pre_cluster_index_[1]) cluster_index_[i] = pre_cluster_index_[0];
		if (cluster_index_[i] > pre_cluster_index_[1]) cluster_index_[i] -= 1;
	}
	cluster_node_count_[pre_cluster_index_[0]] += cluster_node_count_[pre_cluster_index_[1]];
	for (int i = pre_cluster_index_[1]; i < cluster_center_.size() - 1; ++i){
		cluster_center_[i] = cluster_center_[i + 1];
		cluster_node_count_[i] = cluster_node_count_[i + 1];
	}
	cluster_center_.resize(cluster_center_.size() - 1);
	current_cluster_num_ = cluster_center_.size();
	cluster_node_count_.resize(cluster_node_count_.size() - 1);
}

void HierSolver::KmeansCluster() {
	/// Generate random center
	cluster_center_.resize(cluster_num_);
	cluster_node_count_.resize(cluster_num_);
	std::vector< bool > is_selected;
	is_selected.resize(value_.size(), false);
	srand((unsigned int)time(0));
	for (int i = 0; i < cluster_num_; ++i) {
		int temp_index;
		do {
			temp_index = (int)((float)rand() / RAND_MAX * (value_.size() - 1) + 0.499);
		} while (is_selected[temp_index]);

		cluster_center_[i] = value_[temp_index];
		is_selected[temp_index] = true;
	}
	for (int i = 0; i < cluster_index_.size(); ++i) cluster_index_[i] = -1;

	/// Update cluster
	bool is_center_updated;
	int iteration_count = 0;
	do {
		is_center_updated = false;
		for (int i = 0; i < cluster_index_.size(); ++i) {
			float min_dis_index = -1;
			float min_dis = 1e20;
			for (int j = 0; j < cluster_num_; ++j) {
				float temp_dis = 0;
				for (int k = 0; k < weight_.size(); ++k)
					temp_dis += abs(cluster_center_[j][k] - value_[i][k]) * weight_[k];
				if (temp_dis < min_dis) {
					min_dis = temp_dis;
					min_dis_index = j;
				}
			}
			if (cluster_index_[i] != min_dis_index) {
				is_center_updated = true;
				cluster_index_[i] = min_dis_index;
			}
		}
		if (!is_center_updated) break;

		/// update cluster center
		for (int i = 0; i < cluster_center_.size(); ++i)
			memset(cluster_center_[i].data(), 0, cluster_center_[i].size() * sizeof(float));
		memset(cluster_node_count_.data(), 0, cluster_node_count_.size() * sizeof(int));
		for (int i = 0; i < cluster_index_.size(); ++i) {
			cluster_node_count_[cluster_index_[i]]++;
			for (int j = 0; j < weight_.size(); ++j)
				cluster_center_[cluster_index_[i]][j] += value_[i][j];
		}
		for (int i = 0; i < cluster_num_; ++i)
			if (cluster_node_count_[i] != 0) {
				for (int j = 0; j < weight_.size(); ++j)
					cluster_center_[i][j] /= cluster_node_count_[i];
			}

		iteration_count++;
	} while (is_center_updated && iteration_count < 3);
}