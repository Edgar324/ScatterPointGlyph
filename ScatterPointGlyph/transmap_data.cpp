#include "transmap_data.h"
#include "scatter_point_dataset.h"
#include <vtkDelaunay2D.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkLine.h>
#include <vtkCellLocator.h>
#include <vtkSmartPointer.h>
#include <vtkDelaunay2D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCellArray.h>
#include <vtkCell.h>


TransMapData::TransMapData() {
	min_point_num = 5;
}

TransMapData::~TransMapData() {
}

void TransMapData::ProcessData() {
	std::vector< CLeaf* > leaf_nodes;
	for (int i = 0; i < cluster_num; ++i) {
		CLeaf* leaf = new CLeaf;
		leaf->set_level(0);
		leaf->center_pos.resize(2, 0);
		leaf->average_values.resize(this->var_num, 0);

		cluster_node_map.insert(std::map< int, CNode* >::value_type(i, leaf));
		leaf_nodes.push_back(leaf);
	}

	cluster_point_count.resize(this->cluster_num);
	cluster_point_count.assign(this->cluster_num, 0);

	for (int i = 0; i < cluster_index.size(); ++i) {
		int temp_index = cluster_index[i];
		if (temp_index < 0) continue;
		leaf_nodes[temp_index]->center_pos[0] += dataset->original_point_pos[i][0];
		leaf_nodes[temp_index]->center_pos[1] += dataset->original_point_pos[i][1];
		for (int j = 0; j < dataset->weights.size(); ++j)
			leaf_nodes[temp_index]->average_values[j] += (dataset->original_point_values[i][j] - dataset->original_value_ranges[j][0]) / (dataset->original_value_ranges[j][1] - dataset->original_value_ranges[j][0]);
		cluster_point_count[temp_index]++;
	}

	for (int i = 0; i < this->cluster_num; ++i)
		if (cluster_point_count[i] != 0) {
			leaf_nodes[i]->center_pos[0] /= cluster_point_count[i];
			leaf_nodes[i]->center_pos[1] /= cluster_point_count[i];
			for (int j = 0; j < dataset->weights.size(); ++j)
				leaf_nodes[i]->average_values[j] /= cluster_point_count[i];
		}
	
	std::vector< std::vector< int > > temp_index;
	temp_index.resize(this->cluster_num);
	for (int i = 0; i < leaf_nodes.size(); ++i) temp_index[i].push_back(i);

	level_one_colors.clear();
	level_zero_colors.clear();

	for (int i = 0; i < leaf_nodes.size(); ++i)
		if (cluster_point_count[i] >= min_point_num) {
			level_one_nodes.push_back(leaf_nodes[i]);
			level_one_colors.push_back(this->cluster_reprentative_color[3 * i]);
			level_one_colors.push_back(this->cluster_reprentative_color[3 * i + 1]);
			level_one_colors.push_back(this->cluster_reprentative_color[3 * i + 2]);
		} else {
			level_zero_nodes.push_back(leaf_nodes[i]);
			level_zero_colors.push_back(this->cluster_reprentative_color[3 * i]);
			level_zero_colors.push_back(this->cluster_reprentative_color[3 * i + 1]);
			level_zero_colors.push_back(this->cluster_reprentative_color[3 * i + 2]);
		}

	this->UpdateConnectingStatus();
}

void TransMapData::UpdateConnectingStatus() {
	vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < level_one_nodes.size(); ++i) {
		points->InsertNextPoint(level_one_nodes[i]->center_pos[0], level_one_nodes[i]->center_pos[1], 0);
	}
	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer< vtkDelaunay2D >::New();
	delaunay->SetInputData(polydata);
	delaunay->SetTolerance(0.00001);
	delaunay->SetBoundingTriangulation(false);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

	this->node_connecting_status.resize(this->level_one_nodes.size());
	for (int k = 0; k < this->level_one_nodes.size(); ++k) 
		this->node_connecting_status[k].resize(this->level_one_nodes.size(), false);
	for (int k = 0; k < triangle_out->GetNumberOfPolys(); ++k){
		vtkCell* cell = triangle_out->GetCell(k);
		int id1 = cell->GetPointId(0);
		int id2 = cell->GetPointId(1);
		int id3 = cell->GetPointId(2);
		this->node_connecting_status[id1][id2] = true;
		this->node_connecting_status[id2][id1] = true;
		this->node_connecting_status[id2][id3] = true;
		this->node_connecting_status[id3][id2] = true;
		this->node_connecting_status[id1][id3] = true;
		this->node_connecting_status[id3][id1] = true;
	}
}