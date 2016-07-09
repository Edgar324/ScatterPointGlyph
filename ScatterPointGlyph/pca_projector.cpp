/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "pca_projector.h"
#include "vtkSmartPointer.h"
#include "vtkLine.h"
#include "vtkLineSource.h"
#include "vtkTransform.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPCAStatistics.h"
#include "vtkStringArray.h"
#include "vtkTable.h"

PcaProjector::PcaProjector() {

}

PcaProjector::~PcaProjector() {

}

void PcaProjector::Evaluate(vector<vector<float>>& vals) {
    vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
    for (int i = 0; i < vals.size(); ++i) {
        vtkSmartPointer<vtkDoubleArray> xArray = vtkSmartPointer<vtkDoubleArray>::New();
        xArray->SetNumberOfComponents(1);
        char buffer[10];
        itoa(i, buffer, 10);
        xArray->SetName(buffer);
        for (int j = 0; j < vals[i].size(); ++j)
            xArray->InsertNextValue(vals[i][j]);
        datasetTable->AddColumn(xArray);
    }

    vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();

    pcaStatistics->SetInputData( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable );

    for (int i = 0; i < vals.size(); ++i) {
        char buffer[10];
        itoa(i, buffer, 10);

        pcaStatistics->SetColumnStatus(buffer, 1);
    }
    pcaStatistics->RequestSelectedColumns();
    pcaStatistics->SetDeriveOption(true);
    pcaStatistics->Update();

    vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
    

    main_axis_.resize(2);
    for (int i = 0; i < main_axis_.size(); ++i) {
        main_axis_[i].resize(vals.size());
        pcaStatistics->GetEigenvector(i, eigenvectors);
        for (int j = 0; j < main_axis_[i].size(); ++j)
            main_axis_[i][j] = eigenvectors->GetValue(j);
    }
}

void PcaProjector::Apply(vector<vector<float>>& vals, vector<vector<float>>& proj_pos) {
    int record_num = vals[0].size();
    int var_num = vals.size();
    proj_pos.resize(2);
    proj_pos[0].resize(record_num);
    proj_pos[1].resize(record_num);
    for (int i = 0; i < record_num; ++i) {
        proj_pos[0][i] = 0;
        for (int j = 0; j < var_num; ++j)
            proj_pos[0][i] += main_axis_[0][j] * vals[j][i];
        proj_pos[1][i] = 0;
        for (int j = 0; j < var_num; ++j)
            proj_pos[1][i] += main_axis_[1][j] * vals[j][i];
    }
}