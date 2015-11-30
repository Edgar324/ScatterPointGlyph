#include "gestalt_candidate_set.h"
#include <iostream>
#include <QTime>
#include "scatter_point_dataset.h"
#include "Scagnostics.h"

GestaltCandidateSet::GestaltCandidateSet(ScatterPointDataset* data)
	: dataset_(data) {
}

GestaltCandidateSet::~GestaltCandidateSet() {

}

void GestaltCandidateSet::InitSiteData(int num) {
	site_num = num;
	site_point_num.resize(site_num);
	site_center_pos.resize(site_num);
	site_average_value.resize(site_num);
	is_site_labeled.resize(site_num);
	for (int i = 0; i < site_num; ++i) {
		site_point_num[i] = 0;
		site_center_pos[i].resize(2);
		site_center_pos[i][0] = 0;
		site_center_pos[i][1] = 0;
		site_average_value[i] = 0;
		is_site_labeled[i] = false;
	}

	point_site_id.resize(dataset_->point_pos.size());
}

void GestaltCandidateSet::ExtractGestaltCandidates(float dis_thresh) {
	gestalt_candidates.resize(site_num);

	std::vector< bool > is_reached;
	is_reached.resize(site_num);

	std::vector< float > node_distance;
	node_distance.resize(site_num);

	for (int i = 0; i < site_num; ++i) {
		node_distance.assign(site_num, 1e10);
		is_reached.assign(site_num, false);
		node_distance[i] = 0;

		for (int j = 0; j < site_num; ++j) {
			float min_dist = 1e20;
			int min_index = -1;
			for (int k = 0; k < site_num; ++k)
				if (!is_reached[k] && node_distance[k] < min_dist && node_distance[k] < dis_thresh) {
					min_dist = node_distance[k];
					min_index = k;
				}
			if (min_index == -1) break;
			is_reached[min_index] = true;

			for (int k = 0; k < site_num; ++k)
				if (!is_reached[k] && site_connecting_status[min_index][k]) {
					float temp_dis = sqrt(pow(site_center_pos[min_index][0] - site_center_pos[k][0], 2) + pow(site_center_pos[min_index][1] - site_center_pos[k][1], 2));
					if (node_distance[k] > node_distance[min_index] + temp_dis) node_distance[k] = node_distance[min_index] + temp_dis;
				}
		}
		
		gestalt_candidates[i].clear();
		gestalt_candidates[i].push_back(i);
		for (int j = 0; j < site_num; ++j)
			if (is_reached[j] && j != i) gestalt_candidates[i].push_back(j);
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

	QTime time = QTime::currentTime();

	gestalt_scagnostics.resize(1);
	gestalt_scagnostics[0].resize(9);

	std::vector< double > x;
	std::vector< double > y;
	x.resize(dataset_->point_pos.size());
	y.resize(dataset_->point_pos.size());
	for (int i = 0; i < dataset_->point_pos.size(); ++i) {
		x[i] = dataset_->point_pos[i][0];
		y[i] = dataset_->point_pos[i][1];
	}

	Binner b;
	BinnedData bdata = b.binHex(dataset_->point_pos.size(), x.data(), y.data(), 20);
	Triangulation dt;

	double* r = dt.compute(bdata, false);
	//gestalt_scagnostics[i][9] = bdata.n;
	memcpy(gestalt_scagnostics[0].data(), r, sizeof(double) * 9);

	int t = QTime::currentTime().msecsTo(time);

	std::cout << "Time cost: " << t << std::endl;
}