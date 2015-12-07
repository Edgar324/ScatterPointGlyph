#ifndef LINKED_TREE_H_
#define LINKED_TREE_H_

#include <vector>
#include "tree_common.h"

class ScatterPointDataset;
class GestaltProcessor2;
class GestaltCandidateSet;


class LinkedTree
{
public:
	LinkedTree(ScatterPointDataset* data);
	~LinkedTree();

	void ConstructOnOctree(int level_num, float thre);
	void ConstructOnKmeans(int level_num, int basic_cnum);
	void ConstructDirectly(int level_num);
	void UpdateConnectingStatus();

	std::vector< std::vector< CNode*> > tree_nodes;
	std::vector< std::vector< bool > > site_connecting_status;

private:
	ScatterPointDataset* dataset_;

	void ClearData();
	void ConstructLevel(int level);
	void VtkTriangulation();
};

#endif
