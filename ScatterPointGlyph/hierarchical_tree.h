#ifndef HIERARCHICAL_TREE_H_
#define HIERARCHICAL_TREE_H_

#include "tree_common.h"

class ScatterPointDataset;

class HierarchicalTree : public TreeCommon
{
public:
	HierarchicalTree(ScatterPointDataset* data);
	~HierarchicalTree();

	enum DistanceType {
		MIN_DISTANCE = 0x0,
		CENTER_DISTANCE,
		MAX_DISTANCE
	};

	void SetExpectedClusterNum(int num);
	void SetDistanceType(DistanceType type);

protected:
	virtual void GenerateClusters();
	virtual void SplitNode(CBranch* node);

private:
	int expected_cluster_num_;
	DistanceType type_;
};

#endif