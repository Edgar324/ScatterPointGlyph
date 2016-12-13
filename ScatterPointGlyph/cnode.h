#ifndef CNODE_H_
#define CNODE_H_

#include <vector>
using namespace std;

class MultivariateDataset;

class CNode
{
public:
	CNode(MultivariateDataset* dataset_t);
	virtual ~CNode();

	enum NodeType {
		LEAF = 0x0,
		BRANCH
	};

	virtual NodeType type() = 0;

	int level() { return level_; }
	void set_level(int l) { level_ = l; }
	int id() { return id_; }
    unsigned int point_count() { return point_count_; }

    const vector<double>& mean_pos() { return mean_pos_; }
    const vector<double>& mean_values() { return mean_values_; }
    const vector<double>& std_deviations() { return std_deviations_; }

    // Update the statistics for the cluster of data
    void Update(bool value = true, bool pos = true, bool std_dev = true);

protected:
    MultivariateDataset* mv_dataset_ = NULL;

    // Statistics for the cluster of data
	unsigned int point_count_ = 0;
	vector<double> mean_pos_;
	vector<double> mean_values_;
	vector<double> std_deviations_;

private:
	int level_;
	int id_;

	static int max_id_;
};

class CLeaf : public CNode
{
public:
	CLeaf(MultivariateDataset* dataset_t, vector<int>& point_ids_t);
    CLeaf(MultivariateDataset* dataset_t, int point_id_t);
	virtual ~CLeaf();

	virtual CNode::NodeType type() { return CNode::LEAF; }

    const vector<int>& point_ids() { return point_ids_; }

protected:
	vector<int> point_ids_;
};

class CBranch : public CNode
{
public:
	CBranch(MultivariateDataset* dataset_t, vector<CNode*>& children_t);
	virtual ~CBranch();

	virtual CNode::NodeType type() { return CNode::BRANCH; }

    const vector<CNode*>& children() { return children_; }
    void ClearChildren();
    void AppendChild(CNode* node);
    void AppendChildren(vector<CNode*>& nodes);

protected:
	vector<CNode*> children_;
};

#endif