/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef VOLUME_DATA_READER_H_
#define VOLUME_DATA_READER_H_

#include "point_data_reader.h"

class VolumeDataReader : public PointDataReader
{
public:
    VolumeDataReader();
    ~VolumeDataReader();

    virtual ScatterPointDataset* LoadFile(const char* file_name); 
};

#endif