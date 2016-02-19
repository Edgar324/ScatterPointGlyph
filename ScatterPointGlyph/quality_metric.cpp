#include "quality_metric.h"
#include <fstream>

#include "scatter_point_dataset.h"

QualityMetric::QualityMetric() {

}

QualityMetric::~QualityMetric() {

}

void QualityMetric::GenerateQualityMeasures(TreeCommon* tree) {
	int max_level = tree->GetMaxLevel();

	quality_measures_.clear();
	quality_measures_.resize(max_level);

	for (int i = 0; i < max_level; ++i)
		this->GenerateLvelMeasure(tree, i, quality_measures_[i]);
}

void QualityMetric::GenerateLvelMeasure(TreeCommon* tree, int level, std::vector< float >& measures) {
	std::vector< CNode* > cluster_nodes;
	tree->GetClusterResult(level, cluster_nodes);

	int cluster_num = cluster_nodes.size();

	ScatterPointDataset* dataset = tree->data();
	float node_rate = (float)cluster_nodes.size() / dataset->point_num;
	float nnm = 0;

	for (int i = 0; i < dataset->point_num; ++i) {
		float temp_dis = 1e10;
		float x = dataset->point_pos[i][0];
		float y = dataset->point_pos[i][1];
		for (int j = 0; j < cluster_num; ++j) {
			float dis = sqrt(pow(x - cluster_nodes[j]->center_pos[0], 2) + pow(y - cluster_nodes[j]->center_pos[1], 2));
			if (dis < temp_dis) temp_dis = dis;
		}
		nnm += temp_dis;
	}
	nnm = 1.0 - nnm / dataset->point_num;

	measures.resize(2);
	measures[0] = node_rate;
	measures[1] = nnm;
}

void QualityMetric::SaveMeasures(const char* file_path) {
	std::ofstream output(file_path);

	if (output.good()) {
		output << quality_measures_.size() << std::endl;
		for (int i = 0; i < quality_measures_.size(); ++i)
			output << quality_measures_[i][0] << " " << quality_measures_[i][1] << std::endl;
		output.close();
	}
}