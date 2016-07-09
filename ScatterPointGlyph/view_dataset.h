/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef VIEW_DATASET_H_
#define VIEW_DATASET_H_

#include <QObject>

class ViewDataset : public QObject
{
    Q_OBJECT

public:
    ViewDataset() {}
    virtual ~ViewDataset() {}

    virtual void Clear() = 0;
    virtual void Modified() {}

signals:
    void DataChanged();
};

#endif