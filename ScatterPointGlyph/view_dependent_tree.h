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

class MultivariateDataset;

class ViewDependentTree : public MultiLabelTree
{
public:
    ViewDependentTree(MultivariateDataset* data);
    virtual ~ViewDependentTree();

    virtual void GetNodes(float left, float right, float bottom, float top, vector<CNode*>& nodes);

    virtual TreeType type() { return TreeCommon::GEO_VIEW_DEPENDENT_MULTI_LABEL_TREE; }
    virtual void ConstructTree(float left, float right, float bottom, float top);

protected:
    float left_, right_, bottom_, top_;

    virtual void BuildLeafs();
    void BuildOnPixels();
    void BuildOnRecords();
};

#endif