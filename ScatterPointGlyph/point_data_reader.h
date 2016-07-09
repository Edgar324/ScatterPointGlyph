/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef POINT_DATA_READER_H_
#define POINT_DATA_READER_H_

class ScatterPointDataset;

class PointDataReader
{
public:
    PointDataReader();
    virtual ~PointDataReader();

    virtual ScatterPointDataset* LoadFile(const char* file_name);

private:
};

#endif