#ifndef NCUT_TREE_H_
#define NCUT_TREE_H_

#include "tree_common.h"

class NCutTree : public TreeCommon
{
public:
	NCutTree(ScatterPointDataset* data);
	~NCutTree();

	void SetUncertaintyThreshold(float un_threshold);

protected:
	float un_threshold_;

	virtual void GenerateClusters();
	virtual void SplitNode(CBranch* node);
};

#endif