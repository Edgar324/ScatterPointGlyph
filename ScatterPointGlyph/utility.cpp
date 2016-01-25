#include "utility.h"
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

#include "cnode.h"

Utility::Utility() {
}

Utility::~Utility() {

}

void Utility::Sort(std::vector<int>& index_one, std::vector<int>& index_two) {
	for (int i = 0; i < index_one.size() - 1; ++i)
		for (int j = i + 1; j < index_one.size(); ++j)
			if (index_one[i] > index_one[j]){
				int temp_value = index_one[j];
				index_one[j] = index_one[i];
				index_one[i] = temp_value;

				int temp_index = index_two[i];
				index_two[i] = index_two[j];
				index_two[j] = temp_index;
			}
}

void Utility::VtkTriangulation(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status, float& min_edge_length) {
	connecting_status.resize(nodes.size());
	for (int i = 0; i < nodes.size(); ++i) {
		connecting_status[i].resize(nodes.size());
		connecting_status[i].assign(nodes.size(), false);
	}

	vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < nodes.size(); ++i) {
		points->InsertNextPoint(nodes[i]->center_pos[0], nodes[i]->center_pos[1], 0);
	}
	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->SetTolerance(0.00001);
	delaunay->SetBoundingTriangulation(false);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

	min_edge_length = 0.0;
	int edge_count = 0;

	vtkIdTypeArray* idarray = triangle_out->GetPolys()->GetData();
	float min_dis = 1e10;
	for (int i = 0; i < triangle_out->GetNumberOfPolys(); ++i){
		vtkCell* cell = triangle_out->GetCell(i);
		int id1 = cell->GetPointId(0);
		int id2 = cell->GetPointId(1);
		int id3 = cell->GetPointId(2);
		connecting_status[id1][id2] = true;
		connecting_status[id2][id1] = true;
		connecting_status[id2][id3] = true;
		connecting_status[id3][id2] = true;
		connecting_status[id1][id3] = true;
		connecting_status[id3][id1] = true;

		float temp_dis;
		temp_dis = sqrt(pow(nodes[id1]->center_pos[0] - nodes[id2]->center_pos[0], 2)
			+ pow(nodes[id1]->center_pos[1] - nodes[id2]->center_pos[1], 2));
		if (temp_dis < min_dis) min_dis = temp_dis;

		temp_dis = sqrt(pow(nodes[id1]->center_pos[0] - nodes[id3]->center_pos[0], 2)
			+ pow(nodes[id1]->center_pos[1] - nodes[id3]->center_pos[1], 2));
		if (temp_dis < min_dis) min_dis = temp_dis;

		temp_dis = sqrt(pow(nodes[id3]->center_pos[0] - nodes[id2]->center_pos[0], 2)
			+ pow(nodes[id3]->center_pos[1] - nodes[id2]->center_pos[1], 2));
		if (temp_dis < min_dis) min_dis = temp_dis;
	}
	min_edge_length = min_dis;

#ifdef DEBUG_ON
	/*for (int i = 0; i < connecting_status.size(); ++i){
	bool is_connecting = false;
	for (int j = 0; j < connecting_status[i].size(); ++j) is_connecting = is_connecting || connecting_status[i][j];
	assert(is_connecting);
	}*/
#endif // DEBUG_ON
}