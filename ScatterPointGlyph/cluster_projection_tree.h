/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef CLUSTER_PROJECTION_TREE_H_
#define CLUSTER_PROJECTION_TREE_H_

#include "multi_label_tree.h"
#include <vector>
using namespace std;

class ScatterPointDataset;

class ClusterProjectionTree : public MultiLabelTree
{
public:
    ClusterProjectionTree(ScatterPointDataset* data);
    virtual ~ClusterProjectionTree();

    virtual TreeType type() { return TreeCommon::CLUSTER_PROJECTION_TREE; }
    virtual void ConstructTree(float left, float right, float bottom, float top, float glyph_radius);

private:
    vector<CNode*> leaf_nodes;
};

#endif