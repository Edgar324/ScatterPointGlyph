#include "cnode.h"

int CNode::max_id_ = 0;

CNode::CNode() {
	this->id_ = max_id_++;
}

CNode::~CNode() {

}

CLeaf::CLeaf() : CNode() {
}

CLeaf::~CLeaf() {

}

CBranch::CBranch() : CNode() {
}

CBranch::~CBranch() {
	for (int i = 0; i <linked_nodes.size(); ++i)
		delete linked_nodes[i];
	linked_nodes.clear();
}