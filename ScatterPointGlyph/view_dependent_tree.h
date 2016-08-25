/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef VIEW_DEPENDENT_TREE_H_
#define VIEW_DEPENDENT_TREE_H_

#include "multi_label_tree.h"
#include <vector>
using namespace std;

class ScatterPointDataset;

class ViewDependentTree : public MultiLabelTree
{
public:
    ViewDependentTree(ScatterPointDataset* data);
    virtual ~ViewDependentTree();

    virtual TreeType type() { return TreeCommon::VIEW_DEPENDENT_TREE; }
    virtual void ConstructTree(float left, float right, float bottom, float top, float glyph_radius);

private:
    vector<CNode*> leaf_nodes;
};

#endif