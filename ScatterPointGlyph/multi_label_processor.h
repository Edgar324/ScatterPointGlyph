#ifndef MULTI_LABEL_PROCESSOR_H_
#define MULTI_LABEL_PROCESSOR_H_

#include <vector>

class MultiLabelProcessor 
{
public:
	MultiLabelProcessor();
	~MultiLabelProcessor();

	void SetData(std::vector<std::vector<float>>& pos, std::vector<std::vector<float>>& value, std::vector < float>& weights, std::vector<std::vector<bool>>& edges);
	void SetLabelEstimationRadius(float radius);
	void SetSampleNumber(int num);
	void GenerateCluster();
	std::vector<int>& GetResultLabel();

private:
	int point_num_;
	std::vector<std::vector<float>> point_pos_;
	std::vector<std::vector<float>> point_value_;
	std::vector<float> var_weights_;
	std::vector<std::vector<bool>> edges_;
	std::vector<std::vector<float>> edge_weights_;
	std::vector<std::vector<float>> point_dis_;

	const float data_dis_scale_ = { 0.0 };

	float max_radius_;
	int sample_num_;

    float average_value_dis_;

	std::vector<int> result_label_;

	std::vector<std::vector<int>> estimated_models_;

	std::vector<std::vector<float>> data_cost_;
	std::vector<std::vector<float>> smooth_cost_;
	std::vector<float> label_cost_;

	void UpdateEnergyCost();
	void ExtractEstimatedModels();
};

#endif