#ifndef MULTI_LABEL_TREE_H_
#define MULTI_LABEL_TREE_H_

#include "tree_common.h"

class MultiLabelProcessor;

class MultiLabelTree : public TreeCommon
{
public:
	MultiLabelTree(ScatterPointDataset* data);
	virtual ~MultiLabelTree();

    virtual void AutoConstructTree(float std_dev_threshold);
    virtual TreeType type() { return TreeCommon::MULTI_LABEL_TREE; }

protected:
	void SplitNode(CBranch* node);
    void SplitOnSlic(CBranch* node);
    void DirectSplit(CBranch* node);

    void SplitPoints(vector<vector<float>>& pos, vector<vector<float>>& value, float radius, vector<int>& clusters);

	MultiLabelProcessor* processor_ = NULL;

	float max_radius_threshold_ = 0.5;
	float factor_ = 2.0;
    int SLIC_PIXEL_THRESHOLD = 500;
    const int EXPECTED_CLUSTER_NUM = 20;
};

#endif