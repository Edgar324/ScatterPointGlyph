/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "glyph_dataset.h"
#include "glyph_object.h"

GlyphDataset::GlyphDataset() {

}

GlyphDataset::~GlyphDataset() {
    this->Clear();
}

void GlyphDataset::AddGlyphObject(GlyphObject* object) {
    if (object == NULL) return;
    map<int, GlyphObject*>::iterator iter = glyph_object_map_.find(object->cluster_id());
    if (iter == glyph_object_map_.end()) {
        glyph_object_map_.insert(map<int, GlyphObject*>::value_type(object->cluster_id(), object));
    }
}

GlyphObject* GlyphDataset::GetGlyphObject(int cluster_id) {
    map<int, GlyphObject*>::iterator iter = glyph_object_map_.find(cluster_id);
    if (iter != glyph_object_map_.end()) 
        return iter->second;
    else
        return NULL;
}

void GlyphDataset::GetAllGlyphObjects(vector<GlyphObject*>& objects) {
    objects.clear();
    map<int, GlyphObject*>::iterator iter = glyph_object_map_.begin();
    while (iter != glyph_object_map_.end()) {
        objects.push_back(iter->second);
        iter++;
    }
}

void GlyphDataset::Clear() {
    map<int, GlyphObject*>::iterator iter = glyph_object_map_.begin();
    while (iter != glyph_object_map_.end()) {
        delete iter->second;
        iter++;
    }
    glyph_object_map_.clear();
}

void GlyphDataset::Modified() {
    emit DataUpdated();
}