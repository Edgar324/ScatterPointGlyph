/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "var_selection_widget.h"
#include <QtCore/QStringList>

VariableSelectionWidget::VariableSelectionWidget() {
    item_model_ = new QStandardItemModel(0, 3, this);
    this->setModel(item_model_);

    connect(this, SIGNAL(clicked(QModelIndex)), this, SLOT(OnItemClicked(QModelIndex)));

    this->UpdateWidget();
}

VariableSelectionWidget::~VariableSelectionWidget() {

}

void VariableSelectionWidget::SetData(vector<QString>& names, vector<QColor>& colors) {
    var_names_ = names;
    var_colors_ = colors;

    this->UpdateWidget();
}

void VariableSelectionWidget::GetSelection(vector<bool>& is_selected) {
    is_selected.clear();
    for (int i = 0; i < item_model_->rowCount(); ++i) {
        if (item_model_->data(item_model_->item(i, 0)->index(), Qt::CheckStateRole) == 0)
            is_selected.push_back(false);
        else
            is_selected.push_back(true);
    }
}

void VariableSelectionWidget::UpdateWidget() {
    item_model_->clear();

    QStringList headers;
    headers << "Visible" << "Name" << "Color";
    item_model_->setHorizontalHeaderLabels(headers);

    DisableColumnEidtDelegate *dced = new DisableColumnEidtDelegate();
    this->setItemDelegateForColumn(0, dced );
    this->setItemDelegateForColumn(2, dced );
    this->setColumnWidth(0, 80);
    this->setColumnWidth(1, 80);
    this->setColumnWidth(2, 80);


    for (int i = 0; i < var_names_.size(); ++i) {
        item_model_->insertRow(i);
        item_model_->setData(item_model_->index(i, 0), Qt::Checked, Qt::CheckStateRole);
        item_model_->setData(item_model_->index(i, 1), var_names_[i], Qt::DisplayRole);
        item_model_->item(i, 1)->setTextAlignment(Qt::AlignCenter);
        item_model_->setData(item_model_->index(i, 2), QVariant(var_colors_[i]), Qt::DecorationRole);
    }
}

void VariableSelectionWidget::OnItemClicked(QModelIndex index) {
    if (index.column() == 0) {
        if (item_model_->data(index, Qt::CheckStateRole) == 0)
            item_model_->setData(index, Qt::Checked, Qt::CheckStateRole);
        else
            item_model_->setData(index, Qt::Unchecked, Qt::CheckStateRole);
        emit SelectionChanged();
    }
}