#include "cnode.h"

int CNode::max_id_ = 0;

CNode::CNode() : level_(-1) {
	this->id_ = max_id_++;
	this->level_ = -1;

	this->radius = 0;
	this->is_expanded = false;
	this->is_highlighted = false;

	this->color = QColor(128, 128, 128, 200);
	this->hstart = 0;
	this->hend = 1.0;

	this->point_count = 0;
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
	for (int i = 0; i < linked_nodes.size(); ++i)
		delete linked_nodes[i];
	linked_nodes.clear();
}