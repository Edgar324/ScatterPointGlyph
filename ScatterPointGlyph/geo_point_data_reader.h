/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University

 Contact: liaohs082@gmail.com

 All rights reserved.
 */

#ifndef GEO_POINT_DATA_READER_H_
#define GEO_POINT_DATA_READER_H_

#include "point_data_reader.h"

class GeoPointDataReader : public PointDataReader 
{
public:
    GeoPointDataReader();
    virtual ~GeoPointDataReader();

    virtual ScatterPointDataset* LoadFile(const char* file_name); 
};

#endif