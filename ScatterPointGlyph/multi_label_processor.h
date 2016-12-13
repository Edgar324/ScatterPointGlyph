#ifndef MULTI_LABEL_PROCESSOR_H_
#define MULTI_LABEL_PROCESSOR_H_

#include <vector>
using namespace std;

class MultiLabelProcessor 
{
public:
	MultiLabelProcessor(double label_cost_rate = 1.0, double data_dis_scale = 0.5);
	~MultiLabelProcessor();

    void GenerateCluster(const vector<vector<double>>& pos, const vector<vector<double>>& value, 
        const vector<double>& weights, const vector<vector<bool>>& edges,
        double radius, vector<vector<int>>& clusters);
    void GenerateCluster(const vector<vector<double>>& value, 
        const vector<double>& weights, const vector<vector<bool>>& edges,
        double radius, vector<vector<int>>& clusters);

private:
    double label_cost_rate_, data_dis_scale_;

    double max_radius_;
    vector<vector<double>> dis_mat_;
    vector<vector<int>> estimated_models_;

    vector<vector<double>> data_cost_;

    float average_value_dis_;
	vector<vector<double>> smooth_cost_;

	vector<double> label_cost_;

    void GenerateCluster(const vector<vector<bool>>& edges, vector<vector<int>>& clusters);
	void ExtractEstimatedModels(const vector<vector<bool>>& edges);
};

#endif