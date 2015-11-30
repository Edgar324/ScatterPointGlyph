#ifndef LINKED_TREE_H_
#define LINKED_TREE_H_

#include <vector>

class ScatterPointDataset;
class GestaltProcessor2;
class GestaltCandidateSet;

class Node 
{
public:
	Node();
	virtual ~Node();

	enum NodeType {
		LEAF = 0x0,
		BRANCH,
		UNKNOWN
	};

	virtual NodeType type() { return type_; }
	int level() { return level_; }
	void set_level(int l) { level_ = l; }

	int seq_index;

protected:
	NodeType type_;
	int level_;
};

class Leaf : public Node
{
public:
	Leaf();
	~Leaf();

	std::vector< int > linked_points;
};

class Branch : public Node
{
public:
	Branch();
	~Branch();

	std::vector< Node* > linked_nodes;
};

class LinkedTree
{
public:
	LinkedTree(ScatterPointDataset* data);
	~LinkedTree();

	void ConstructOnOctree(int level_num, float thre);
	void ConstructOnKmeans(int level_num, int basic_cnum);
	void ConstructDirectly(int level_num);
	void UpdateConnectingStatus();

	std::vector< std::vector< Node*> > tree_nodes;
	std::vector< std::vector< bool > > site_connecting_status;

private:
	ScatterPointDataset* dataset_;

	void ClearData();
	void ConstructLevel(int level);
	void VtkTriangulation();
};

#endif
