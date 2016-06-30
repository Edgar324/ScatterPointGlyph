/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef VAR_SELECTION_WIDGET_H_
#define VAR_SELECTION_WIDGET_H_

#include <QTableView>
#include <QStyledItemDelegate>
#include <QtGui/QColor>
#include <QStandardItemModel>
#include <vector>
using namespace std;

class DisableColumnEidtDelegate : public QStyledItemDelegate   
{   
public:   
    DisableColumnEidtDelegate(QObject *parent = 0): QStyledItemDelegate(parent) { }   
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option,   
        const QModelIndex &index) const   
    {
        return NULL;   
    }
};

class VarSelectionWidget : public QTableView
{
    Q_OBJECT

public:
    VarSelectionWidget();
    ~VarSelectionWidget();

    void SetData(vector<QString>& names, vector<QColor>& colors);
    void GetSelection(vector<bool>& is_selected);

signals:
    void SelectionChanged();

private:
    vector<QString> var_names_;
    vector<QColor> var_colors_;

    QStandardItemModel* item_model_;

    void UpdateWidget();

    private slots:
    void OnItemClicked(QModelIndex);
};

#endif