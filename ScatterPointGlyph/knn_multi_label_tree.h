/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef KNN_MULTI_LABEL_TREE_H_
#define KNN_MULTI_LABEL_TREE_H_

#include "tree_common.h"

class MultiLabelProcessor;
class CBranch;
class MultivariateDataset;

class KnnMultiLabelTree : public TreeCommon
{
public:
	KnnMultiLabelTree(MultivariateDataset* data, int nearest_neighbors = 5);
	virtual ~KnnMultiLabelTree();

    virtual TreeType type() { return TreeCommon::KNN_MULTI_LABEL_TREE; }

    virtual void GetNodes(int level, vector<CNode*>& level_nodes);
    virtual void GetNodes(float left, float right, float bottom, float top, vector<CNode*>& nodes);

    virtual void SplitNode(CBranch* node);
    void SplitNodes(int level, vector<CBranch*>& nodes);

    virtual void AutoConstructTree(float std_dev_threshold);
    virtual void ConstructTree(float left, float right, float bottom, float top);
    virtual void Clear();

protected:
    MultiLabelProcessor* processor_ = NULL;

    vector<vector<CNode*>> exploration_path_;

    vector<int> current_records_;

    int nearest_neighbors_ = 5;
    double radius_scale_ = 0.7;

    const int SCALE_NODE_SIZE = 2000;
    const int LEAF_RESO_SIZE = 1000;
    const int SLIC_RESO_SIZE = 200;

    const int SLIC_PIXEL_THRESHOLD = 500;
    const int EXPECTED_CLUSTER_NUM = 40;

    virtual void BuildLeafs();
    void BuildOnKwayPartitions();
    void BuildOnRecords();
};

#endif