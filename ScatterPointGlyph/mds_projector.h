/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef MDS_PROJECTOR_H_
#define MDS_PROJECTOR_H_

#include <vector>
using namespace std;

class MdsProjector
{
public:
    MdsProjector();
    ~MdsProjector();

    // Execute MDS on the normalized multi-variate data
	// The implementation of MDS is from Quan Wang.
	// More information about this function can be retrieved from "SimpleMatrix.h" or
	// https://sites.google.com/site/simpmatrix/
    static bool Project(vector<vector<float>>& values, vector<vector<float>>& proj_values);
};

#endif