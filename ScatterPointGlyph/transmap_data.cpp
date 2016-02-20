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
	min_point_num = 3;
}

TransMapData::~TransMapData() {
}

void TransMapData::ClearData() {
	cluster_num = 0;
	var_num = 0;
	cluster_node_map.clear();

	level_one_nodes.clear();
}

void TransMapData::ProcessData() {
	this->cluster_node_map.clear();

	for (int i = 0; i < cluster_nodes.size(); ++i)
		this->cluster_node_map.insert(std::map< int, CNode* >::value_type(i, cluster_nodes[i]));

	for (int i = 0; i < cluster_nodes.size(); ++i)
		if (cluster_nodes[i]->point_count >= min_point_num) {
			level_one_nodes.push_back(cluster_nodes[i]);
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
	for (int k = 0; k < this->level_one_nodes.size(); ++k){
		this->node_connecting_status[k].resize(this->level_one_nodes.size(), false);
		this->node_connecting_status[k].assign(this->level_one_nodes.size(), false);
	}
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