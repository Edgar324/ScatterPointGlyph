/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#ifndef GLYPH_DESIGN_DIALOG_H_
#define GLYPH_DESIGN_DIALOG_H_

#include <QtWidgets/QDialog>
#include <QtWidgets/QCheckBox>
#include <QtCore/QString>
#include <vector>
#include "ui_glyph_design_dialog.h"

class QVTKWidget;
class vtkRenderer;
class vtkActor;
class vtkPolyDataMapper;
class vtkPolyData;
class vtkTextActor3D;

class GlyphDesignDialog : public QDialog
{
    Q_OBJECT
public:
    GlyphDesignDialog();
    ~GlyphDesignDialog();

    void SetVariables(std::vector< QString >& names);
    void GetSelectedVariables(std::vector< int >& index);

private:
    Ui::GlyphDesignDialog ui_;
    QVTKWidget* vtkwidget_;
    vtkRenderer* main_renderer_;

    std::vector< QCheckBox* > variable_checkbox_;
    std::vector< float > mean_values_;
    std::vector< float > msd_values_;

    vtkActor* glyph_actor_;
    vtkPolyData* glyph_poly_;
    vtkPolyDataMapper* glyph_mapper_;
    std::vector< vtkTextActor3D* > text_actors_;

    void UpdateGlyph();

    private slots:
    void OnVariableSelectionChanged();
};

#endif
