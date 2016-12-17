#include "cnode.h"
#include <queue>
#include <iostream>
using namespace std;
#include "multivariate_dataset.h"

int CNode::max_id_ = 0;

CNode::CNode(MultivariateDataset* dataset_t) 
    : mv_dataset_(dataset_t), id_(max_id_++) {
}

CNode::~CNode() {

}

void CNode::GetNodeRecords(vector<int>& ids) {
    ids.clear();

	std::queue<CNode*> node_queue;
	node_queue.push(this);

	while (node_queue.size() != 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast< CBranch* >(temp_node);
			if (branch != NULL) {
                const vector<CNode*>& children = branch->children();
				for (int i = 0; i < children.size(); ++i)
					node_queue.push(children[i]);
			}
		}
		else if (temp_node->type() == CNode::LEAF) {
			CLeaf* leaf = dynamic_cast< CLeaf* >(temp_node);
			if (leaf != NULL) {
                const vector<int>& record_ids = leaf->record_ids();
				for (int i = 0; i < record_ids.size(); ++i)
					ids.push_back(record_ids[i]);
			}
		}
	}
}

void CNode::Update(bool update_value/* = true*/, bool update_pos/* = true*/) {

	std::vector<int> record_ids;
	this->GetNodeRecords(record_ids);

    if (record_ids.size() == 0) {
        cout << "Error: Node " << this->id() << " have 0 points!" << endl;
        return;
    }

    int var_num = mv_dataset_->var_num();

    this->point_count_ = record_ids.size();
    this->mean_values_.resize(var_num, 0);
	this->std_deviations_.resize(var_num, 0);
	this->mean_pos_.resize(2, 0);

    if (update_value) 
        mv_dataset_->GetAverageValues(record_ids, mean_values_, std_deviations_);

    if (update_pos) {
        double x, y;
        mv_dataset_->GetAveragePos(record_ids, x, y);
        mean_pos_[0] = x;
        mean_pos_[1] = y;
    }
}

CLeaf::CLeaf(MultivariateDataset* dataset_t, vector<int>& point_ids_t) 
    : CNode(dataset_t), record_ids_(point_ids_t) {
}

CLeaf::CLeaf(MultivariateDataset* dataset_t, int point_id_t) 
    : CLeaf(dataset_t, vector<int>(1, point_id_t)) {
    
}

CLeaf::~CLeaf() {

}

void CLeaf::ClearRecords() {
    record_ids_.clear();
}

void CLeaf::AppendRecord(int record_id) {
    record_ids_.push_back(record_id);
}

CBranch::CBranch(MultivariateDataset* dataset_t, vector<CNode*>& children_t)
    : CNode(dataset_t), children_(children_t) {
    for (int i = 0; i < children_.size(); i++)
        children_[i]->set_level(this->level() + 1);
}

CBranch::CBranch(MultivariateDataset* dataset_t, vector<CLeaf*>& children_t) 
    : CNode(dataset_t) {
    children_.resize(children_t.size());
    for (int i = 0; i < children_.size(); i++) {
        children_[i] = children_t[i];
        children_[i]->set_level(this->level() + 1);
    }
}

CBranch::CBranch(MultivariateDataset* dataset_t) 
    : CNode(dataset_t) {

}

CBranch::~CBranch() {
	for (int i = 0; i < children_.size(); ++i) 
		delete children_[i];
	children_.clear();
}

void CBranch::ClearChildren() {
    children_.clear();
}

void CBranch::AppendChild(CNode* node) {
    children_.push_back(node);
}

void CBranch::AppendChildren(vector<CNode*>& nodes) {
    int csize = children_.size();
    children_.resize(csize + nodes.size());
    for (int i = 0; i < nodes.size(); i++)
        children_[csize + i] = nodes[i];
}
