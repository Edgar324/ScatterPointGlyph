#ifndef CNODE_H_
#define CNODE_H_

#include <vector>
#include <iostream>
using namespace std;
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

    // Bounding rect of each node
    float left = 0, right = 0, bottom = 0, top = 0;
    vector<float> center_pos;

    // Statistics of each node
	int point_count = 0;
	vector<float> mean_pos;
	vector<float> mean_values;
	vector<float> std_deviations;

	// The parameter for generating this node using multi-label clustering method
	float radius = 0;

    //
    float average_dis = 0;

    bool is_expanded = false;
    bool is_visible = false;

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

	vector<int> linked_points;
};

class CBranch : public CNode
{
public:
	CBranch();
	virtual ~CBranch();

	virtual CNode::NodeType type() { return CNode::BRANCH; }

	vector<CNode*> linked_nodes;
};

#endif