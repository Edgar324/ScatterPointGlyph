#include "gestalt_candidate_set.h"
#include <iostream>
#include <QTime>
#include "scatter_point_dataset.h"
#include "Scagnostics.h"

GestaltCandidateSet::GestaltCandidateSet() {
}

GestaltCandidateSet::~GestaltCandidateSet() {

}

void GestaltCandidateSet::ExtractGestaltCandidates(float dis_thresh) {
	int cluster_num = clusters.size();

	gestalt_candidates.resize(cluster_num);
	gestalt_cluster_index.resize(cluster_num);

	std::vector< bool > is_reached;
	is_reached.resize(cluster_num);

	std::vector< float > node_distance;
	node_distance.resize(cluster_num);

	for (int i = 0; i < cluster_num; ++i) {
		node_distance.assign(cluster_num, 1e10);
		is_reached.assign(cluster_num, false);
		node_distance[i] = 0;

		for (int j = 0; j < cluster_num; ++j) {
			float min_dist = 1e20;
			int min_index = -1;
			for (int k = 0; k < cluster_num; ++k)
				if (!is_reached[k] && node_distance[k] < min_dist && node_distance[k] < dis_thresh) {
					min_dist = node_distance[k];
					min_index = k;
				}
			if (min_index == -1) break;
			is_reached[min_index] = true;

			for (int k = 0; k < cluster_num; ++k)
				if (!is_reached[k] && cluster_connecting_status[min_index][k]) {
					float temp_dis = sqrt(pow(clusters[min_index]->center_pos[0] - clusters[k]->center_pos[0], 2) + pow(clusters[min_index]->center_pos[1] - clusters[k]->center_pos[1], 2));
					if (node_distance[k] > node_distance[min_index] + temp_dis) node_distance[k] = node_distance[min_index] + temp_dis;
				}
		}
		
		gestalt_candidates[i].clear();
		gestalt_cluster_index[i].clear();
		for (int j = 0; j < cluster_num; ++j)
			if (is_reached[j]) {
				gestalt_cluster_index[i].push_back(j);
				for (int k = 0; k < site_nodes.size(); ++k)
					if (basic_node_index[k] == j) gestalt_candidates[i].push_back(k);
			}
	}
}

void GestaltCandidateSet::InitSiteData() {
	this->cluster_connecting_status.resize(this->clusters.size());
	for (int i = 0; i < this->clusters.size(); ++i) {
		this->cluster_connecting_status[i].resize(this->clusters.size());
		this->cluster_connecting_status[i].assign(this->clusters.size(), false);
	}

	for (int i = 0; i < site_nodes.size() - 1; ++i)
		for (int j = i + 1; j < site_nodes.size(); ++j)
			if (site_connecting_status[i][j]) {
				int cluster_one = basic_node_index[i];
				int cluster_two = basic_node_index[j];
				cluster_connecting_status[cluster_one][cluster_two] = true;
			}
}

void GestaltCandidateSet::ExtractScagnostics() {
	//gestalt_scagnostics.resize(gestalt_candidates.size());
	//for (int i = 0; i < gestalt_candidates.size(); ++i) {
	//	std::vector< double > x;
	//	std::vector< double > y;
	//	x.resize(gestalt_candidates[i].size());
	//	y.resize(gestalt_candidates[i].size());

	//	for (int j = 0; j < gestalt_candidates[i].size(); ++j) {
	//		int pIndex = gestalt_candidates[i][j];
	//		x[j] = dataset_->point_pos[pIndex][0];
	//		y[j] = dataset_->point_pos[pIndex][1];
	//	}

	//	gestalt_scagnostics[i].resize(9);

	//	Binner b;
	//	BinnedData bdata = b.binHex(gestalt_candidates[i].size(), x.data(), y.data(), 20);
	//	Triangulation dt;

	//	double* r = dt.compute(bdata, false);
	//	//gestalt_scagnostics[i][9] = bdata.n;
	//	memcpy(gestalt_scagnostics[i].data(), r, sizeof(double) * 9);
	//	/*for (int i = 0; i < bdata.n; i++) {
	//		results[10 + 0 * bdata.n + i] = bdata.x[i];
	//		results[10 + 1 * bdata.n + i] = bdata.y[i];
	//		results[10 + 2 * bdata.n + i] = bdata.counts[i];
	//	}*/
	//}

	//QTime time = QTime::currentTime();

	//gestalt_scagnostics.resize(1);
	//gestalt_scagnostics[0].resize(9);

	//std::vector< double > x;
	//std::vector< double > y;
	//x.resize(dataset_->point_pos.size());
	//y.resize(dataset_->point_pos.size());
	//for (int i = 0; i < dataset_->point_pos.size(); ++i) {
	//	x[i] = dataset_->point_pos[i][0];
	//	y[i] = dataset_->point_pos[i][1];
	//}

	//Binner b;
	//BinnedData bdata = b.binHex(dataset_->point_pos.size(), x.data(), y.data(), 20);
	//Triangulation dt;

	//double* r = dt.compute(bdata, false);
	////gestalt_scagnostics[i][9] = bdata.n;
	//memcpy(gestalt_scagnostics[0].data(), r, sizeof(double) * 9);

	//int t = QTime::currentTime().msecsTo(time);

	//std::cout << "Time cost: " << t << std::endl;
}