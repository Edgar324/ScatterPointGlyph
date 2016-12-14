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

    const int SCALE_NODE_SIZE = 2000;
    const int LEAF_RESO_SIZE = 1000;
    const int SLIC_RESO_SIZE = 200;

    const int SLIC_PIXEL_THRESHOLD = 500;
    const int EXPECTED_CLUSTER_NUM = 40;

    bool is_pixel_used_ = false;

    virtual void BuildLeafs();
    void BuildOnPixels();
    void BuildOnRecords();

    void SplitOnPixels(const vector<CNode*>& children, vector<vector<int>>& clusters);
    void SplitOnRecords(const vector<CNode*>& children, vector<vector<int>>& clusters);
    void SplitPoints(vector<vector<double>>& pos, vector<vector<double>>& value, double radius, vector<vector<int>>& clusters);
};

#endif