#include "tree_common.h"
#include <queue>
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

TreeCommon::TreeCommon(ScatterPointDataset* data)
	: dataset_(data),
	average_edge_length_(0) {

	root_ = new CBranch;
}

TreeCommon::~TreeCommon() {

}

void TreeCommon::GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index) {

}

void TreeCommon::run() {

}

void TreeCommon::ConstructOnOctree(float thre) {

}

void TreeCommon::ConstructOnKmeans(int basic_cnum) {

}

void TreeCommon::ConstructDirectly() {
	leaf_nodes_.resize(dataset_->point_pos.size());
	for (int i = 0; i < dataset_->point_pos.size(); ++i) {
		CLeaf* temp_leaf = new CLeaf();
		temp_leaf->center_pos = dataset_->point_pos[i];
		temp_leaf->average_values = dataset_->point_values[i];
		temp_leaf->linked_points.push_back(i);
		leaf_nodes_[i] = temp_leaf;
		leaf_nodes_[i]->seq_index = i;
		leaf_nodes_[i]->set_level(0);
	}

	if (dataset_->is_structured_data) {
		node_connecting_status_.resize(dataset_->point_pos.size());
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			node_connecting_status_[i].resize(dataset_->point_pos.size());
			node_connecting_status_[i].assign(node_connecting_status_[i].size(), false);
		}
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			int x = i % dataset_->w;
			int y = i / dataset_->w;
			if (x == dataset_->w - 1 || y == dataset_->h - 1) continue;
			int pre_one = y * dataset_->w + x + 1;
			int pre_two = (y + 1) * dataset_->w + x;
			node_connecting_status_[pre_one][i] = true;
			node_connecting_status_[i][pre_one] = true;
			node_connecting_status_[pre_two][i] = true;
			node_connecting_status_[i][pre_two] = true;
		}
	}
	else {
		VtkTriangulation();
	}
}

void TreeCommon::GenerateCluster(int min_pixel_radius) {

}

void TreeCommon::Traverse(int level, std::vector< CNode* >& nodes) {
	nodes.clear();

	std::queue< CNode* > node_queue;
	node_queue.push(root_);
	while (node_queue.size() != 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->level() == level) {
			nodes.push_back(temp_node);
		}
		else if (temp_node->level() > level) {
			if (temp_node->type() == CNode::BRANCH) {
				CBranch* branch = dynamic_cast<CBranch*>(temp_node);
				if (branch != NULL) {
					for (int i = 0; i < branch->linked_nodes.size(); ++i)
						node_queue.push(branch->linked_nodes[i]);
				}
			}
		}
	}
}

void TreeCommon::Traverse(CNode* node, std::vector< int >& linked_points) {
	std::queue< CNode* > node_queue;
	node_queue.push(node);

	while (node_queue.size() != 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast< CBranch* >(temp_node);
			if (branch != NULL) {
				for (int i = 0; i < branch->linked_nodes.size(); ++i)
					node_queue.push(branch->linked_nodes[i]);
			}
		}
		else if (temp_node->type() == CNode::LEAF) {
			CLeaf* leaf = dynamic_cast< CLeaf* >(temp_node);
			if (leaf != NULL) {
				for (int i = 0; i < leaf->linked_points.size(); ++i)
					linked_points.push_back(leaf->linked_points[i]);
			}
		}
	}
}

void TreeCommon::VtkTriangulation() {
	vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < leaf_nodes_.size(); ++i) {
		CLeaf* leaf_node = dynamic_cast<CLeaf*>(leaf_nodes_[i]);
		points->InsertNextPoint(leaf_node->center_pos[0], leaf_node->center_pos[1], 0);
	}
	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->SetTolerance(1e-30);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

	int site_num = leaf_nodes_.size();
	node_connecting_status_.resize(site_num);
	for (int i = 0; i < site_num; ++i) node_connecting_status_[i].resize(site_num, false);

	average_edge_length_ = 0.0;

	vtkIdTypeArray* idarray = triangle_out->GetPolys()->GetData();
	for (int i = 0; i < triangle_out->GetNumberOfPolys(); ++i){
		double* temp = idarray->GetTuple(i * 3);
		int id1 = (int)(temp[0]);
		temp = idarray->GetTuple(i * 3 + 1);
		int id2 = (int)(temp[0]);
		temp = idarray->GetTuple(i * 3 + 2);
		int id3 = (int)(temp[0]);
		node_connecting_status_[id1][id2] = true;
		node_connecting_status_[id2][id1] = true;
		node_connecting_status_[id2][id3] = true;
		node_connecting_status_[id3][id2] = true;
		node_connecting_status_[id1][id3] = true;
		node_connecting_status_[id3][id1] = true;

		average_edge_length_ += sqrt(pow(leaf_nodes_[id1]->center_pos[0] - leaf_nodes_[id2]->center_pos[0], 2)
			+ pow(leaf_nodes_[id1]->center_pos[1] - leaf_nodes_[id2]->center_pos[1], 2));
		average_edge_length_ += sqrt(pow(leaf_nodes_[id1]->center_pos[0] - leaf_nodes_[id3]->center_pos[0], 2)
			+ pow(leaf_nodes_[id1]->center_pos[1] - leaf_nodes_[id3]->center_pos[1], 2));
		average_edge_length_ += sqrt(pow(leaf_nodes_[id3]->center_pos[0] - leaf_nodes_[id2]->center_pos[0], 2)
			+ pow(leaf_nodes_[id3]->center_pos[1] - leaf_nodes_[id2]->center_pos[1], 2));
	}
	average_edge_length_ /= (3 * triangle_out->GetNumberOfPolys());
	// For temperary usage
	average_edge_length_ = 0.05;

	for (int i = 0; i < node_connecting_status_.size(); ++i){
		bool is_connecting = false;
		for (int j = 0; j < node_connecting_status_[i].size(); ++j) is_connecting = is_connecting || node_connecting_status_[i][j];
		assert(is_connecting);
	}
}