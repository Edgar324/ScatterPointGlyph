#ifndef QUALITY_METRIC_H_
#define QUALITY_METRIC_H_

#include "tree_common.h"
#include <vector>

class QualityMetric
{
public:
	QualityMetric();
	~QualityMetric();

	void GenerateQualityMeasures(TreeCommon* tree);
	void SaveMeasures(const char* file_path);

    void GenerateLvelMeasure(TreeCommon* tree, int level, std::vector<float>& measures);

private:
	std::vector<std::vector<float>> quality_measures_;
};

#endif