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
#include <vtkCellArray.h>
#include <vtkCell.h>

#include "scatter_point_dataset.h"

int CNode::max_id_ = 0;

CNode::CNode() : type_(UNKNOWN), level_(-1) {
	this->id = max_id_++;
	is_expanded = true;
}

CNode::~CNode() {

}

CLeaf::CLeaf() : CNode() {
	type_ = CNode::LEAF;
}

CLeaf::~CLeaf() {

}

CBranch::CBranch() : CNode() {
	type_ = CNode::BRANCH;
}

CBranch::~CBranch() {

}

CNode* CBranch::FindNearestNode(float x, float y) {
	return NULL;
}

CNode* CBranch::FindNearestValue(std::vector< float >& values) {
	return NULL;
}

TreeCommon::TreeCommon(ScatterPointDataset* data)
	: dataset_(data),
	min_edge_length_(0),
	max_level_(0) {

	root_ = new CBranch;
}

TreeCommon::~TreeCommon() {

}

void TreeCommon::GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index) {

}

void TreeCommon::GetClusterResult(float dis_per_pixel, int& cluster_num, std::vector< int >& cluster_index) {

}

void TreeCommon::run() {

}

void TreeCommon::ConstructOnOctree(float thre) {

}

void TreeCommon::ConstructOnKmeans(int basic_cnum) {

}

void TreeCommon::ConstructOnRandomSample(int sample_num) {
	srand((unsigned int)time(0));
	std::vector< bool > is_used;
	is_used.resize(dataset_->point_pos.size(), false);

	root_->linked_nodes.clear();

	float temp_min = 1e10;
	for (int i = 0; i < dataset_->point_pos.size(); ++i)
		if (!is_used[i]) {
			is_used[i] = true;
			CLeaf* temp_leaf = new CLeaf();
			std::vector< int > neighbour_list;
			neighbour_list.push_back(i);
			for (int j = i + 1; j < dataset_->point_pos.size(); ++j){
				float dis = sqrt(pow(dataset_->point_pos[j][0] - dataset_->point_pos[i][0], 2) + pow(dataset_->point_pos[j][1] - dataset_->point_pos[i][1], 2));
				if (dis < 1e-5) {
					neighbour_list.push_back(j);
					is_used[j] = true;
				}
				if (dis < temp_min) temp_min = dis;
			}
			temp_leaf->linked_points = neighbour_list;

			temp_leaf->center_pos.resize(2);
			temp_leaf->center_pos[0] = 0;
			temp_leaf->center_pos[1] = 0;
			temp_leaf->average_values.resize(dataset_->weights.size());
			for (int j = 0; j < dataset_->weights.size(); ++j)
				temp_leaf->average_values[j] = 0;
			for (int j = 0; j < neighbour_list.size(); ++j) {
				temp_leaf->center_pos[0] += dataset_->point_pos[neighbour_list[j]][0];
				temp_leaf->center_pos[1] += dataset_->point_pos[neighbour_list[j]][1];
				for (int k = 0; k < dataset_->weights.size(); ++k)
					temp_leaf->average_values[k] += dataset_->point_values[neighbour_list[j]][k];
			}
			temp_leaf->center_pos[0] /= neighbour_list.size();
			temp_leaf->center_pos[1] /= neighbour_list.size();
			for (int j = 0; j < dataset_->weights.size(); ++j)
				temp_leaf->average_values[j] /= neighbour_list.size();

			root_->linked_nodes.push_back(temp_leaf);
			temp_leaf->set_level(0);
		}
}

void TreeCommon::ConstructDirectly() {
	float min_dis_threshold = 1e-5;
	std::vector< bool > is_used;
	is_used.resize(dataset_->point_pos.size(), false);

	root_->linked_nodes.clear();

	float temp_min = 1e10;
	for (int i = 0; i < dataset_->point_pos.size(); ++i) 
		if (!is_used[i]) {
			is_used[i] = true;
			CLeaf* temp_leaf = new CLeaf();
			std::vector< int > neighbour_list;
			neighbour_list.push_back(i);
			for (int j = i + 1; j < dataset_->point_pos.size(); ++j){
				float dis = sqrt(pow(dataset_->point_pos[j][0] - dataset_->point_pos[i][0], 2) + pow(dataset_->point_pos[j][1] - dataset_->point_pos[i][1], 2));
				if (dis < 1e-5) {
					neighbour_list.push_back(j);
					is_used[j] = true;
				}
				if (dis < temp_min) temp_min = dis;
			}
			temp_leaf->linked_points = neighbour_list;

			temp_leaf->center_pos.resize(2);
			temp_leaf->center_pos[0] = 0;
			temp_leaf->center_pos[1] = 0;
			temp_leaf->average_values.resize(dataset_->weights.size());
			for (int j = 0; j < dataset_->weights.size(); ++j)
				temp_leaf->average_values[j] = 0;
			for (int j = 0; j < neighbour_list.size(); ++j) {
				temp_leaf->center_pos[0] += dataset_->point_pos[neighbour_list[j]][0];
				temp_leaf->center_pos[1] += dataset_->point_pos[neighbour_list[j]][1];
				for (int k = 0; k < dataset_->weights.size(); ++k)
					temp_leaf->average_values[k] += dataset_->point_values[neighbour_list[j]][k];
			}
			temp_leaf->center_pos[0] /= neighbour_list.size();
			temp_leaf->center_pos[1] /= neighbour_list.size();
			for (int j = 0; j < dataset_->weights.size(); ++j)
				temp_leaf->average_values[j] /= neighbour_list.size();
			temp_leaf->point_count = temp_leaf->linked_points.size();

			root_->linked_nodes.push_back(temp_leaf);
			temp_leaf->set_level(0);
		}

	min_edge_length_ = 1e10;
	for (int i = 0; i < root_->linked_nodes.size() - 1; ++i)
		for (int j = i + 1; j < root_->linked_nodes.size(); ++j) {
			float temp_dis;
			temp_dis = sqrt(pow(root_->linked_nodes[i]->center_pos[0] - root_->linked_nodes[j]->center_pos[0], 2)
				+ pow(root_->linked_nodes[i]->center_pos[1] - root_->linked_nodes[j]->center_pos[1], 2));
			if (temp_dis < min_edge_length_) min_edge_length_ = temp_dis;
		}
	/*if (dataset_->is_structured_data) {
		node_connecting_status_.resize(dataset_->point_pos.size());
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			node_connecting_status_[i].resize(dataset_->point_pos.size());
			node_connecting_status_[i].assign(node_connecting_status_[i].size(), false);
		}
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			int x = dataset_->sample_index[i] % dataset_->w;
			int y = dataset_->sample_index[i] / dataset_->w;
			if (x == dataset_->w - 1 || y == dataset_->h - 1) continue;
			int pre_one = -1;
			if (dataset_->node_sample_map.find(y * dataset_->w + x + 1) != dataset_->node_sample_map.end())
				pre_one = dataset_->node_sample_map[y * dataset_->w + x + 1];
			if (pre_one == -1) continue;
			int pre_two = -1;
			if (dataset_->node_sample_map.find((y + 1) * dataset_->w + x) != dataset_->node_sample_map.end())
				pre_two = dataset_->node_sample_map[(y + 1) * dataset_->w + x];
			if (pre_two == -1) continue;
			node_connecting_status_[pre_one][i] = true;
			node_connecting_status_[i][pre_one] = true;
			node_connecting_status_[pre_two][i] = true;
			node_connecting_status_[i][pre_two] = true;
		}

		this->min_edge_length_ = sqrt(pow(leaf_nodes_[0]->center_pos[0] - leaf_nodes_[1]->center_pos[0], 2)
			+ pow(leaf_nodes_[0]->center_pos[1] - leaf_nodes_[1]->center_pos[1], 2));
	}
	else {
		
	}*/
}

void TreeCommon::GenerateCluster(CBranch* node /* = NULL */) {

}

void TreeCommon::Traverse(int level, std::vector< CNode* >& nodes) {
	nodes.clear();

	std::queue< CNode* > node_queue;
	node_queue.push(root_);
	while (node_queue.size() != 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();
		temp_node->is_expanded = true;

		if (temp_node->level() == level) {
			nodes.push_back(temp_node);
			temp_node->is_expanded = false;
		} else {
			if (temp_node->type() == CNode::BRANCH && temp_node->level() < level) {
				CBranch* branch = dynamic_cast<CBranch*>(temp_node);
				bool is_all_leaf = true;
				for (int i = 0; i < branch->linked_nodes.size(); ++i)
					if (branch->linked_nodes[i]->type() != CNode::LEAF) {
						is_all_leaf = false;
						break;
					}
				if (is_all_leaf && level != max_level_) {
					nodes.push_back(temp_node);
					temp_node->is_expanded = false;
				} else {
					for (int i = 0; i < branch->linked_nodes.size(); ++i)
						node_queue.push(branch->linked_nodes[i]);
				}
			}
		}
	}
}

void TreeCommon::Traverse(float radius, std::vector< CNode* >& nodes) {
	nodes.clear();

	std::queue< CNode* > node_queue;
	node_queue.push(root_);
	while (node_queue.size() != 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->radius > radius) {
			nodes.push_back(temp_node);
		}
		else {
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

void TreeCommon::VtkTriangulation(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status) {
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

	min_edge_length_ = 0.0;
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
	min_edge_length_ = min_dis;

#ifdef DEBUG_ON
	for (int i = 0; i < connecting_status.size(); ++i){
		bool is_connecting = false;
		for (int j = 0; j < connecting_status[i].size(); ++j) is_connecting = is_connecting || connecting_status[i][j];
		assert(is_connecting);
	}
#endif // DEBUG_ON
}