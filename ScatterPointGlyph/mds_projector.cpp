/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "mds_projector.h"
#include <assert.h>
#include "SimpleMatrix.h"

MdsProjector::MdsProjector() {

}

MdsProjector::~MdsProjector() {

}

bool MdsProjector::Project(vector<vector<float>>& values, vector<vector<float>>& proj_values) {
    int dim_num = values.size();
    int record_num = values[0].size();
    smat::Matrix<double> *X0 = new smat::Matrix<double>(record_num, dim_num, 0.0);

	for (int i = 0; i < dim_num; ++i) {
		float min_value = 1e30;
		float max_value = -1e30;
		for (int j = 0; j < record_num; ++j) {
			if (values[i][j] > max_value) max_value = values[i][j];
			if (values[i][j] < min_value) min_value = values[i][j];
		}

		assert(max_value - min_value != 0);

		for (int j = 0; j < record_num; ++j) {
			float scale_value = (values[i][j] - min_value) / (max_value - min_value);
			X0->set(j, i, scale_value);
		}
	}

	smat::Matrix<double> *D = new smat::Matrix<double>(record_num, record_num, 0.0);
	for (int i = 0; i < record_num - 1; ++i){
		D->set(i, i, 0);
		for (int j = i + 1; j < record_num; ++j) {
			float dis = 0;
			for (int k = 0; k < dim_num; ++k)
				dis += pow((X0->get(i, k) - X0->get(j, k)), 2) + 1e-20;
			dis = sqrt(dis);
			assert(dis != 0);
			D->set(i, j, dis);
			D->set(j, i, dis);
		}
	}

	int projected_dim = 2;
	int iteration = 40;

	smat::Matrix<double> * X1 = MDS_SMACOF(D, NULL, projected_dim, iteration); // without initialization

	proj_values.resize(2);
    proj_values[0].resize(record_num);
    proj_values[1].resize(record_num);
	for (int i = 0; i < record_num; ++i) {
		proj_values[0][i] = X1->get(i, 0) * 10000;
		proj_values[1][i] = X1->get(i, 1) * 10000;
	}

    return false;
}