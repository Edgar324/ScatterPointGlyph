#include "linked_tree.h"
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
#include "scatter_point_dataset.h"

Node::Node() : type_(UNKNOWN), level_(-1), seq_index(-1) {

}

Node::~Node() {

}

Leaf::Leaf() : Node() {
	type_ = Node::LEAF;
}

Leaf::~Leaf() {

}

Branch::Branch() : Node() {
	type_ = Node::BRANCH;
}

Branch::~Branch() {

}

LinkedTree::LinkedTree(ScatterPointDataset* data) 
	: dataset_(data) {

}

LinkedTree::~LinkedTree() {

}

void LinkedTree::ConstructLevel(int level) {
	// construct data for gestalt candidate set
	
}

void LinkedTree::ConstructOnOctree(int level_num, float thre) {

}

void LinkedTree::ConstructOnKmeans(int level_num, int basic_cnum) {

}

void LinkedTree::ConstructDirectly(int level_num) {
	ClearData();

	tree_nodes.resize(level_num);
	tree_nodes[0].resize(dataset_->original_point_pos.size());
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i) {
		Leaf* temp_leaf = new Leaf();
		temp_leaf->linked_points.push_back(i);
		tree_nodes[0][i] = temp_leaf;
		tree_nodes[0][i]->seq_index = i;
	}

	VtkTriangulation();
}

void LinkedTree::VtkTriangulation() {
	vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < tree_nodes[0].size(); ++i) {
		float x = 0, y = 0;
		Leaf* leaf_node = dynamic_cast<Leaf*>(tree_nodes[0][i]);
		for (int j = 0; j < leaf_node->linked_points.size(); ++j) {
			int point_index = leaf_node->linked_points[j];
			x += dataset_->point_pos[point_index][0];
			y += dataset_->point_pos[point_index][1];
		}
		x /= leaf_node->linked_points.size();
		y /= leaf_node->linked_points.size();

		points->InsertNextPoint(x, y, 0);
	}
	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

	int site_num = tree_nodes[0].size();
	site_connecting_status.resize(site_num);
	for (int i = 0; i < site_num; ++i) site_connecting_status[i].resize(site_num, false);

	vtkIdTypeArray* idarray = triangle_out->GetPolys()->GetData();
	for (int i = 0; i < triangle_out->GetNumberOfPolys(); ++i){
		double* temp = idarray->GetTuple(i * 3);
		int id1 = (int)(temp[0]);
		temp = idarray->GetTuple(i * 3 + 1);
		int id2 = (int)(temp[0]);
		temp = idarray->GetTuple(i * 3 + 2);
		int id3 = (int)(temp[0]);
		site_connecting_status[id1][id2] = true;
		site_connecting_status[id2][id1] = true;
		site_connecting_status[id2][id3] = true;
		site_connecting_status[id3][id2] = true;
		site_connecting_status[id1][id3] = true;
		site_connecting_status[id3][id1] = true;
	}

	for (int i = 0; i < site_connecting_status.size(); ++i){
		bool is_connecting = false;
		for (int j = 0; j < site_connecting_status[i].size(); ++j) is_connecting = is_connecting || site_connecting_status[i][j];
		assert(is_connecting);
	}
}

void LinkedTree::UpdateConnectingStatus() {
	int top_level = tree_nodes.size() - 1;
	while (tree_nodes[top_level].size() == 0) top_level--;
	if (site_connecting_status.size() == tree_nodes[top_level].size()) return;

	std::vector< std::vector< bool > > temp_status;
	temp_status.resize(tree_nodes[top_level].size());
	for (int i = 0; i < tree_nodes[top_level].size(); ++i) {
		temp_status[i].resize(tree_nodes[top_level].size(), false);
	}
	for (int i = 0; i < tree_nodes[top_level].size() - 1; ++i)
		for (int j = i + 1; j < tree_nodes[top_level].size(); ++j) {
			Branch* node_one = dynamic_cast< Branch* >(tree_nodes[top_level][i]);
			Branch* node_two = dynamic_cast< Branch* >(tree_nodes[top_level][j]);
			for (int si = 0; si < node_one->linked_nodes.size(); ++si)
				for (int sj = 0; sj < node_two->linked_nodes.size(); ++sj) {
					temp_status[i][j] = temp_status[i][j] || site_connecting_status[node_one->linked_nodes[si]->seq_index][node_two->linked_nodes[sj]->seq_index];
				}
			temp_status[j][i] = temp_status[i][j];
		}
	site_connecting_status = temp_status;
}

void LinkedTree::ClearData() {

}