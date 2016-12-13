/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef LEAF_BUILDER_H_
#define LEAF_BUILDER_H_

#include <vector>
using namespace std;

class MultivariateDataset;
class CLeaf;

class LeafBuilder
{
public:
    LeafBuilder();
    ~LeafBuilder();

    void Build(MultivariateDataset* dataset, vector<CLeaf*>& leaf_nodes);
};

#endif