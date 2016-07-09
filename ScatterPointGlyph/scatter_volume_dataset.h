/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef SCATTER_VOLUME_DATASET_H_
#define SCATTER_VOLUME_DATASET_H_

#include "scatter_point_dataset.h"

class ScatterVolumeDataset : public ScatterPointDataset
{
public:
    ScatterVolumeDataset(int x, int y, int z);
    ~ScatterVolumeDataset();

    int* dims() { return dims_; }

    vector<int> voxel_index;

    virtual DataType type() { return VOLUME_DATA; }

private:
    int dims_[3];
};

#endif