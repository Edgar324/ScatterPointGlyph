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

class ScatterPointDataset;

class ViewDependentTree : public MultiLabelTree
{
public:
    ViewDependentTree(ScatterPointDataset* data);
    virtual ~ViewDependentTree();

    virtual void AutoConstructTree(float std_dev_threshold);
};

#endif