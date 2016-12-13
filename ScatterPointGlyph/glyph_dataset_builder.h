/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef GLYPH_DATASET_BUILDER_H_
#define GLYPH_DATASET_BUILDER_H_

#include <vector>
using namespace std;

class TreeCommon;
class GlyphObject;
class GlyphDataset;
class MultivariateDataset;
class CNode;

class GlyphDatasetBuilder
{
public:
    GlyphDatasetBuilder() {}
    ~GlyphDatasetBuilder() {}

    static void Build(MultivariateDataset* mv_dataset, vector<int>& selected_var_index,
        TreeCommon* tree, float left, float right, float bottom, float top, float glyph_half_size,
        GlyphDataset* glyph_dataset);

    // Saliency based on the paper:
    // J. Harel, C. Koch, and P. Perona. Graph-based visual saliency. In Advances
    // in neural information processing systems, pages 545¨C552, 2006.
    static void EvaluateSaliency(MultivariateDataset* mv_dataset, vector<CNode*>& nodes, vector<float>& node_saliency);
};

#endif