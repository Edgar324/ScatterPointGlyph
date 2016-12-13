#include "cnode.h"
#include "multivariate_dataset.h"

int CNode::max_id_ = 0;

CNode::CNode(MultivariateDataset* dataset_t) 
    : mv_dataset_(dataset_t), id_(max_id_++) {
}

CNode::~CNode() {

}

void CNode::Update(bool value/* = true*/, bool pos/* = true*/, bool std_dev/* = true*/) {
    /// TODO: Update the statistics
}

CLeaf::CLeaf(MultivariateDataset* dataset_t, vector<int>& point_ids_t) 
    : CNode(dataset_t), point_ids_(point_ids_t) {
}

CLeaf::CLeaf(MultivariateDataset* dataset_t, int point_id_t) 
    : CLeaf(dataset_t, vector<int>(1, point_id_t)) {
    
}

CLeaf::~CLeaf() {

}

CBranch::CBranch(MultivariateDataset* dataset_t, vector<CNode*>& children_t) 
    : CNode(dataset_t), children_(children_t) {
    for (int i = 0; i < children_.size(); i++)
        children_[i]->set_level(this->level() + 1);
}

CBranch::~CBranch() {
	for (int i = 0; i < children_.size(); ++i)
		delete children_[i];
	children_.clear();
}

void CBranch::ClearChildren() {

}

void CBranch::AppendChild(CNode* node) {

}

void CBranch::AppendChildren(vector<CNode*>& nodes) {

}
