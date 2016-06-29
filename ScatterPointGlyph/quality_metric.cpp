#include "quality_metric.h"
#include <fstream>

#include "scatter_point_dataset.h"

QualityMetric::QualityMetric() {

}

QualityMetric::~QualityMetric() {

}

void QualityMetric::GenerateQualityMeasures(TreeCommon* tree) {
	int max_level = tree->GetMaxLevel();

    if (max_level > 6) max_level = 6;

	quality_measures_.clear();
	quality_measures_.resize(max_level);

	for (int i = 0; i < max_level; ++i)
		this->GenerateLvelMeasure(tree, i, quality_measures_[i]);
}

void QualityMetric::GenerateLvelMeasure(TreeCommon* tree, int level, std::vector<float>& measures) {
	std::vector<CNode*> cluster_nodes;
	tree->GetClusterResult(level, cluster_nodes);

	int cluster_num = cluster_nodes.size();

	// nearest neighbor evaluation
	ScatterPointDataset* dataset = tree->data();
	float node_rate = (float)cluster_nodes.size() / dataset->point_num;
	float nnm = 0;

	for (int i = 0; i < dataset->point_num; ++i) {
		float temp_dis = 1e10;
		int temp_index = -1;
		float x = dataset->normalized_point_pos[i][0];
		float y = dataset->normalized_point_pos[i][1];
		for (int j = 0; j < cluster_num; ++j) {
			float dis = sqrt(pow(x - cluster_nodes[j]->center_pos[0], 2) + pow(y - cluster_nodes[j]->center_pos[1], 2));
			if (dis < temp_dis) {
				temp_dis = dis;
				temp_index = j;
			}
		}
		if (temp_index != -1) {
			float value_dis = 0;
			for (int j = 0; j < dataset->var_num; ++j)
				value_dis += abs(cluster_nodes[temp_index]->average_values[j] - dataset->normalized_point_values[i][j]) * dataset->var_weights[j];
			nnm += value_dis;
		}
	}
	nnm = 1.0 - nnm / dataset->point_num;

	// square evaluation
	float sqrerr = 0.0;
	for (int i = 0; i < cluster_num; ++i) {
		float value_dis = 0;
		if (cluster_nodes[i]->type() != CNode::LEAF) {
			for (int j = 0; j < dataset->var_num; ++j)
				value_dis += cluster_nodes[i]->variable_variances[j] * dataset->var_weights[j];
			sqrerr += value_dis * cluster_nodes[i]->point_count;
		}
	}

	measures.resize(3);
	measures[0] = cluster_nodes.size();
	//measures[1] = nnm;
	measures[2] = sqrerr / tree->root()->point_count;
}

void QualityMetric::SaveMeasures(const char* file_path) {
	std::ofstream output(file_path);

	if (output.good()) {
		output << quality_measures_.size() << std::endl;
		for (int i = 0; i < quality_measures_.size(); ++i)
			output << quality_measures_[i][0] << " " << quality_measures_[i][1] << " " << quality_measures_[i][2] << std::endl;
		output.close();
	}
}