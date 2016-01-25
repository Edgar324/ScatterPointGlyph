#ifndef NCUT_TREE_H_
#define NCUT_TREE_H_

#include "tree_common.h"

class NCutTree : public TreeCommon
{
public:
	NCutTree(ScatterPointDataset* data);
	~NCutTree();

	void SetExpectedClusterNum(int num);
	void SetUncertaintyThreshold(float un_threshold);

protected:
	float un_threshold_;
	float data_dis_scale_;
	int expected_cluster_num_;

	virtual void run();
	virtual void GenerateCluster(CBranch* node);
};

#endif