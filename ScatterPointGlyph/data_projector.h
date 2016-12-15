/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef DATA_PROJECTOR_H_
#define DATA_PROJECTOR_H_

#define EIGEN_RUNTIME_NO_MALLOC

#include <vector>
using namespace std;
#include <Eigen/Dense>
#include "data_core.h"

class DataProjector
{
public:
    DataProjector(ProjectionMethod method_t = tSNE);
    ~DataProjector();

    bool Project(Eigen::MatrixXd& records, Eigen::MatrixXd& projected_pos);
};

#endif