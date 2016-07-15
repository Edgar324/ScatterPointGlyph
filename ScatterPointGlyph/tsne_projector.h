/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef TSNE_PROJECTOR_H_
#define TSNE_PROJECTOR_H_

#include <vector>
using namespace std;

class TsneProjector
{
public:
    TsneProjector();
    ~TsneProjector();

    static void Project(vector<vector<float>>& values, vector<vector<float>>& proj_values);
};

#endif