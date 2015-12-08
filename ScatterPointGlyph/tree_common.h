#ifndef TREE_COMMON_H_
#define TREE_COMMON_H_

#include <vector>
#include <iostream>
#include <QThread>

class ScatterPointDataset;

class CNode
{
public:
	CNode();
	~CNode();

	enum NodeType {
		LEAF = 0x0,
		BRANCH,
		UNKNOWN
	};

	virtual NodeType type() { return type_; }

	std::vector< float > center_pos;
	std::vector< float > average_values;

	int level() { return level_; }
	void set_level(int l) { level_ = l; }
	int seq_index;

protected:
	NodeType type_;
	int level_;
};

class CLeaf : public CNode
{
public:
	CLeaf();
	~CLeaf();

	std::vector< int > linked_points;
};

class CBranch : public CNode
{
public:
	CBranch();
	~CBranch();

	std::vector< CNode* > linked_nodes;
};

class TreeCommon : public QThread
{
public:
	TreeCommon(ScatterPointDataset* data);
	virtual ~TreeCommon();

	void ConstructOnOctree(float thre);
	void ConstructOnKmeans(int basic_cnum);
	void ConstructDirectly();

	virtual void GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index);

protected:
	ScatterPointDataset* dataset_;

	CBranch* root_;
	std::vector< CLeaf* > leaf_nodes_;
	std::vector< std::vector< bool > > node_connecting_status_;

	float average_edge_length_;

	void VtkTriangulation();

	virtual void run();
	virtual void GenerateCluster(int min_pixel_radius = 1);
	void Traverse(CNode* node, std::vector< int >& linked_points);
	void Traverse(int level, std::vector< CNode* >& nodes);
};

#endif