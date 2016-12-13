#ifndef MULTI_LABEL_TREE_H_
#define MULTI_LABEL_TREE_H_

#include "tree_common.h"

class MultiLabelProcessor;

class MultiLabelTree : public TreeCommon
{
public:
	MultiLabelTree(MultivariateDataset* data);
	virtual ~MultiLabelTree();

    virtual TreeType type() { return TreeCommon::GEO_MULTI_LABEL_TREE; }
    virtual void SplitNode(CBranch* node);
    virtual void AutoConstructTree(float std_dev_threshold);
    virtual void ConstructTree(float left, float right, float bottom, float top, float glyph_radius);
    virtual void Clear();

protected:
    MultiLabelProcessor* processor_ = NULL;

	float max_radius_threshold_ = 0.5;
	float factor_ = 2.0;
    int SLIC_PIXEL_THRESHOLD = 500;
    const int EXPECTED_CLUSTER_NUM = 40;

    virtual void BuildLeafs();

    /*void SplitOnSlic(CBranch* node);
    void DirectSplit(CBranch* node);

    void SplitPoints(vector<vector<float>>& pos, vector<vector<float>>& value, float radius, vector<int>& clusters);*/
};

#endif