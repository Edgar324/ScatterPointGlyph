/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "data_projector.h"
#include "tapkee/tapkee.hpp"
using namespace tapkee;

DataProjector::DataProjector(ProjectionMethod method_t/* = tSNE*/) {

}

DataProjector::~DataProjector() {

}

bool DataProjector::Project(Eigen::MatrixXd& records, Eigen::MatrixXd& projected_pos) {
    TapkeeOutput output = tapkee::initialize() 
        .withParameters((method=MultidimensionalScaling, target_dimension=2))
        .embedUsing(records);

    projected_pos = output.embedding;

    return true;
}
