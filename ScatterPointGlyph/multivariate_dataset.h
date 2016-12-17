/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef MULTIVARIATE_DATASET_H_
#define MULTIVARIATE_DATASET_H_

#define EIGEN_RUNTIME_NO_MALLOC

#include <vector>
using namespace std;
#include <Eigen/Dense>
#include <QString>
#include <QColor>
#include "data_core.h"

class DataProjector;

class MultivariateDataset
{
public:
    MultivariateDataset(const char* file_name);
    ~MultivariateDataset();

    enum DataType {
        REGULAR_DATA = 0x0,
        GEO_SPATIAL_DATA,
        IMAGE_DATA,
        VOLUME_DATA
    };

    virtual DataType type() { return REGULAR_DATA; }

    bool is_good() { return is_good_; }

    unsigned int record_num() { return record_num_; }
    unsigned int var_num() { return var_num_; }
    const vector<QString>& var_names() { return var_names_; }
    const vector<double>& var_weights() { return var_weights_; }
    const vector<QColor>& var_colors() { return var_colors_; }
    bool is_projected() { return is_projected_; }

    double* GetRecordPos(int i);
    void GetPosRanges(double& left, double& right, double& bottom, double& top);
    void GetAveragePos(const vector<int>& record_ids, double& x, double& y);
    void GetAverageValues(const vector<int>& record_ids, vector<double>& values, vector<double>& std_deviations);

    // Return the distance of two data records
    // @param record_one record_two The index of two records
    // @param ratio The ratio of data distance in the result, dis = ratio * data_dis + (1.0 - ratio) * geo_dis.
    //              The default ratio is 1.0 for regular multivariate data and can be adjusted for geospatial multivariate data
    virtual double Distance(int record_one, int record_two);
    virtual double Distance(int record_one, int record_two, double ratio);
    
    virtual bool ApplyProjection(DataProjector* projector);

protected:
    // Dataset information
    unsigned int record_num_ = 0;
	unsigned int var_num_ = 0;
	vector<QString> var_names_;
    vector<QColor> var_colors_;

    bool is_projected_ = false;

    // Currently not used, all assigned to 1.0 / var_num_
    vector<double> var_weights_;

    Eigen::MatrixXd records_;
    Eigen::MatrixXd var_ranges_;

    Eigen::MatrixXd projected_pos_;
    Eigen::MatrixXd pos_ranges_;

    double temp_pos_[2];

private:
    // File reading status
    bool is_good_ = false;

    virtual bool LoadFromFile(const char* file_name);
};

#endif