/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef GLYPH_DATASET_H_
#define GLYPH_DATASET_H_

#include <vector>
#include <map>
using namespace std;
#include "view_dataset.h"

class GlyphObject;

class GlyphDataset : public ViewDataset
{
    Q_OBJECT

public:
    GlyphDataset();
    ~GlyphDataset();

    void AddGlyphObject(GlyphObject* object);
    GlyphObject* GetGlyphObject(int cluster_id);
    void GetAllGlyphObjects(vector<GlyphObject*>& objects);

    virtual void Clear();
    virtual void Modified();

signals:
    void DataUpdated();

private:
    map<int, GlyphObject*> glyph_object_map_;
    int highlight_index_ = -1;
};

#endif