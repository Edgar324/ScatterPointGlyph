/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "geo_multivariate_dataset.h"

GeoMultivariateDataset::GeoMultivariateDataset(const char* file_name)
    : MultivariateDataset(file_name) {

}

GeoMultivariateDataset::~GeoMultivariateDataset() {

}

double GeoMultivariateDataset::Distance(int record_one, int record_two, double ratio /*= 0.5*/) {
    return 0.0;
}

bool GeoMultivariateDataset::ApplyProjection(ProjectionMethod method /*= tSNE*/) {
    return false;
}

bool GeoMultivariateDataset::LoadFromFile(const char* file_name) {
    return false;
}

