/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef PARALLEL_DATASET_BUILDER_H_
#define PARALLEL_DATASET_BUILDER_H_

#include <vector>
using namespace std;

class ScatterPointDataset;
class ParallelDataset;

class ParallelDatasetBuilder
{
public:
    ParallelDatasetBuilder() {}
    ~ParallelDatasetBuilder() {}

    static void Build(ScatterPointDataset* point_dataset, vector<int>& selected_var_index, vector<vector<int>>& point_ids, ParallelDataset* parallel_dataset);
};

#endif