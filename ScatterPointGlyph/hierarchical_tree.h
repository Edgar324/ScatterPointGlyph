#ifndef HIERARCHICAL_TREE_H_
#define HIERARCHICAL_TREE_H_

#include "tree_common.h"

class MultivariateDataset;

class HierarchicalTree : public TreeCommon
{
public:
	HierarchicalTree(MultivariateDataset* data);
	~HierarchicalTree();

	enum DistanceType {
		MIN_DISTANCE = 0x0,
		CENTER_DISTANCE,
		MAX_DISTANCE
	};

	void SetExpectedClusterNum(int num);
	void SetDistanceType(DistanceType type);

    virtual void AutoConstructTree(float std_dev_threshold);
    virtual TreeType type() { return TreeCommon::GEO_HIERARCHICAL_TREE; }

protected:
    virtual void SplitNode(CBranch* node);

private:
	int expected_cluster_num_;
	DistanceType type_;
};

#endif