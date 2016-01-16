#ifndef NCUT_TREE_H_
#define NCUT_TREE_H_

#include "tree_common.h"

class NCutTree : public TreeCommon
{
public:
	NCutTree(ScatterPointDataset* data);
	~NCutTree();

	void SetUncertaintyThreshold(float un_threshold);

	virtual void GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index);
	virtual void GetClusterResult(float dis_per_piexl, int& cluster_num, std::vector< int >& cluster_index);
	virtual void GetClusterResult(float radius, std::vector< CNode* >& level_nodes);

protected:
	float un_threshold_;

	virtual void run();
	virtual void GenerateCluster(CBranch* node);
};

#endif