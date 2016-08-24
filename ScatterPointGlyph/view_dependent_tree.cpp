/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "view_dependent_tree.h"

ViewDependentTree::ViewDependentTree(ScatterPointDataset* data)
    : MultiLabelTree(data) {

}

ViewDependentTree::~ViewDependentTree() {

}

void ViewDependentTree::AutoConstructTree(float std_dev_threshold) {

}