/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef GEO_VARIATE_DATASET_H_
#define GEO_VARIATE_DATASET_H_

#include "multivariate_dataset.h"

class GeoMultivariateDataset : public MultivariateDataset
{
public:
    GeoMultivariateDataset(const char* file_name);
    ~GeoMultivariateDataset();

    virtual DataType type() { return GEO_SPATIAL_DATA; }
    virtual double Distance(int record_one, int record_two, double ratio = 0.5);
    virtual bool ApplyProjection(ProjectionMethod method = tSNE);

protected:
    Eigen::MatrixXd record_pos_;
    Eigen::MatrixXd pos_ranges_;

    virtual bool LoadFromFile(const char* file_name);
};

#endif