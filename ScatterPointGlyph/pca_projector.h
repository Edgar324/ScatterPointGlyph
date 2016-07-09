/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef PCA_PROJECTOR_H_
#define PCA_PROJECTOR_H_

#include <vector>
using namespace std;

class PcaProjector
{
public:
    PcaProjector();
    ~PcaProjector();

    void Evaluate(vector<vector<float>>& vals);
    void Apply(vector<vector<float>>& vals, vector<vector<float>>& proj_pos);

private:
    vector<vector<float>> main_axis_;
};

#endif