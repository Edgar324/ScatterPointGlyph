#ifndef TREE_COMMON_H_
#define TREE_COMMON_H_

#include <vector>
#include <iostream>
#include <QThread>

class ScatterPointDataset;
class CBranch;

class CNode
{
public:
	CNode();
	virtual ~CNode();

	enum NodeType {
		LEAF = 0x0,
		BRANCH,
		UNKNOWN
	};

	virtual NodeType type() { return type_; }

	std::vector< float > center_pos;
	std::vector< float > average_values;
	std::vector< float > value_variance;
	std::vector< float > boxplot_upper_bound;
	std::vector< float > boxplot_lower_bound;
	float radius;
	int id;
	int point_count;
	bool is_expanded;
	bool is_highlighted;

	CBranch* parent = NULL;

	int level() { return level_; }
	void set_level(int l) { level_ = l; }

protected:
	NodeType type_;
	int level_;

	static int max_id_;
};

class CLeaf : public CNode
{
public:
	CLeaf();
	virtual ~CLeaf();

	std::vector< int > linked_points;
};

class CBranch : public CNode
{
public:
	CBranch();
	virtual ~CBranch();

	std::vector< CNode* > linked_nodes;

	CNode* FindNearestNode(float x, float y);
	CNode* FindNearestValue(std::vector< float >& values);
};

class TreeCommon : public QThread
{
public:
	TreeCommon(ScatterPointDataset* data);
	virtual ~TreeCommon();

	void ConstructOnOctree(float thre);
	void ConstructOnKmeans(int basic_cnum);
	void ConstructOnRandomSample(int sample_num);
	void ConstructDirectly();

	virtual void GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index);
	virtual void GetClusterResult(float dis_per_pixel, int& cluster_num, std::vector< int >& cluster_index);

	CBranch* root() { return root_; }

protected:
	ScatterPointDataset* dataset_;

	CBranch* root_;
	float min_edge_length_;
	int max_level_;

	virtual void run();
	virtual void GenerateCluster(CBranch* node = NULL);

	void VtkTriangulation(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status);
	void Traverse(CNode* node, std::vector< int >& linked_points);
	void Traverse(int level, std::vector< CNode* >& nodes);
	void Traverse(float radius , std::vector< CNode* >& nodes);
	void ActiveTraverse(std::vector< CNode* >& nodes);
};

#endif