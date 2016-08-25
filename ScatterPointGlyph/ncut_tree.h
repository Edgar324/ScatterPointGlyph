#ifndef NCUT_TREE_H_
#define NCUT_TREE_H_

#include "tree_common.h"

class NCutTree : public TreeCommon
{
public:
	NCutTree(ScatterPointDataset* data);
	~NCutTree();

	virtual void AutoConstructTree(float std_dev_threshold);
    virtual TreeType type() { return TreeCommon::NCUTS_TREE; }

protected:
	virtual void SplitNode(CBranch* node);
};

#endif