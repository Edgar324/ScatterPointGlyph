#ifndef IMMEDIATE_GESTALT_TREE_H_
#define IMMEDIATE_GESTALT_TREE_H_

#include "immediate_tree.h"

class ImmediateGestaltTree : public ImmediateTree
{
public:
	ImmediateGestaltTree(ScatterPointDataset* dataset);
	~ImmediateGestaltTree();

protected:
	virtual void run();
	virtual void GenerateCluster();

};

#endif