/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef SCATTER_GRID_DATASET_H_
#define SCATTER_GRID_DATASET_H_

#include "scatter_point_dataset.h"

class ScatterGridDataset : public ScatterPointDataset
{
public:
    ScatterGridDataset() {
        w = -1;
        h = -1;
    }
    ~ScatterGridDataset() {}

    int w, h;

    virtual DataType type() { return ScatterPointDataset::GRID_DATA; }
};

#endif