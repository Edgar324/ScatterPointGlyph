#ifndef CNODE_H_
#define CNODE_H_

#include <vector>
#include <iostream>
#include <QtGui/QColor>

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

	virtual NodeType type() { return UNKNOWN; }

	CBranch* parent = NULL;

	int point_count;
	std::vector< float > center_pos;
	std::vector< float > average_values;
	std::vector< float > value_variance;

	// The parameter for generating this node using multi-label clustering method
	float radius;

	// Parameters for rendering effect in the main views
	bool is_expanded;
	bool is_highlighted;

	// Representative color for the cluster
	// [hstart, hend] is the hue range for child nodes
	QColor color;
	float hstart, hend;

	int level() { return level_; }
	void set_level(int l) { level_ = l; }
	int id() { return id_; }

private:
	int level_;
	int id_;

	static int max_id_;
};

class CLeaf : public CNode
{
public:
	CLeaf();
	virtual ~CLeaf();

	virtual CNode::NodeType type() { return CNode::LEAF; }

	std::vector< int > linked_points;
};

class CBranch : public CNode
{
public:
	CBranch();
	virtual ~CBranch();

	virtual CNode::NodeType type() { return CNode::BRANCH; }

	std::vector< CNode* > linked_nodes;
	// The placement positions for the linked_nodes
	// This can be generated either by manually or by TSP algorithm
	std::vector< int > sorting_index;
};

#endif