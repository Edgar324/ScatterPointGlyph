/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "scatter_volume_dataset.h"

ScatterVolumeDataset::ScatterVolumeDataset(int x, int y, int z)
    : ScatterPointDataset() {
    this->dims_[0] = x;
    this->dims_[1] = y;
    this->dims_[2] = z;
}

ScatterVolumeDataset::~ScatterVolumeDataset() {

}