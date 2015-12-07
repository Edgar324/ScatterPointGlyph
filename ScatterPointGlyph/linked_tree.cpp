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

CNode::CNode() : type_(UNKNOWN), level_(-1), seq_index(-1) {

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
	std::vector< std::vector< float > > cluster_center;
	std::vector< int > cluster_node_count;
	std::vector< int > cluster_index;
	float dis_scale = 0.5;

	if (dataset_->is_structured_data) {
		site_connecting_status.resize(dataset_->point_pos.size());
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			site_connecting_status[i].resize(dataset_->point_pos.size());
			site_connecting_status[i].assign(site_connecting_status[i].size(), false);
		}
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			int x = i % dataset_->w;
			int y = i / dataset_->w;
			if (x == dataset_->w - 1 || y == dataset_->h - 1) continue;
			int pre_one = y * dataset_->w + x + 1;
			int pre_two = (y + 1) * dataset_->w + x;
			site_connecting_status[pre_one][i] = true;
			site_connecting_status[i][pre_one] = true;
			site_connecting_status[pre_two][i] = true;
			site_connecting_status[i][pre_two] = true;
		}
	}
	else {
		VtkTriangulation();
	}

	/// Generate random center
	cluster_center.resize(basic_cnum);
	cluster_node_count.resize(basic_cnum);

	std::vector< bool > is_selected;
	is_selected.resize(dataset_->point_values.size(), false);

	/// remove duplicated data
	for (int i = 0; i < dataset_->point_pos.size() - 1; ++i)
		if (!is_selected[i]) {
			for (int j = i + 1; j < dataset_->point_pos.size(); ++j) {
				float temp_dis = abs(dataset_->point_pos[i][0] - dataset_->point_pos[j][0]) + abs(dataset_->point_pos[i][1] - dataset_->point_pos[j][1]);
				if (temp_dis <= 1e-20) is_selected[j] = true;
			}
		}

	srand((unsigned int)time(0));
	for (int i = 0; i < basic_cnum; ++i) {
		int temp_index;
		do {
			temp_index = (int)((float)rand() / RAND_MAX * (dataset_->point_values.size() - 1) + 0.499);
		} while (is_selected[temp_index]);

		cluster_center[i].resize(3);
		cluster_center[i][0] = dataset_->point_pos[temp_index][0];
		cluster_center[i][1] = dataset_->point_pos[temp_index][1];
		cluster_center[i][2] = dataset_->point_values[temp_index];
		is_selected[temp_index] = true;
	}

	cluster_index.resize(dataset_->point_pos.size());
	for (int i = 0; i < cluster_index.size(); ++i) cluster_index[i] = -1;

	/// Update cluster
	bool is_center_updated;
	int iteration_count = 0;
	do {
		is_center_updated = false;
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			float min_dis_index = -1;
			float min_dis = 1e20;
			for (int j = 0; j < basic_cnum; ++j) {
				float temp_dis = 0;
				temp_dis += sqrt(pow(dataset_->point_pos[i][0] - cluster_center[j][0], 2) + pow(dataset_->point_pos[i][1] - cluster_center[j][1], 2)) * dis_scale;
				temp_dis += abs(dataset_->point_values[j] - dataset_->point_values[i]) * (1.0 - dis_scale);
				if (temp_dis < min_dis) {
					min_dis = temp_dis;
					min_dis_index = j;
				}
			}
			if (cluster_index[i] != min_dis_index) {
				is_center_updated = true;
				cluster_index[i] = min_dis_index;
			}
		}
		if (!is_center_updated) break;

		/// update cluster center
		for (int i = 0; i < basic_cnum; ++i)
			memset(cluster_center[i].data(), 0, cluster_center[i].size() * sizeof(float));
		memset(cluster_node_count.data(), 0, cluster_node_count.size() * sizeof(int));
		for (int i = 0; i < cluster_index.size(); ++i) {
			cluster_node_count[cluster_index[i]]++;
			cluster_center[cluster_index[i]][0] += dataset_->point_pos[i][0];
			cluster_center[cluster_index[i]][1] += dataset_->point_pos[i][1];
			cluster_center[cluster_index[i]][2] += dataset_->point_values[i];
		}
		for (int i = 0; i < basic_cnum; ++i)
			if (cluster_node_count[i] != 0) {
				for (int j = 0; j < cluster_center[i].size(); ++j)
					cluster_center[i][j] /= cluster_node_count[i];
			}

		iteration_count++;
	} while (is_center_updated && iteration_count < 30);

	for (int i = 0; i < basic_cnum; ++i) {
		if (cluster_node_count[i] == 0) {
			while (cluster_node_count[i] == 0 && i < basic_cnum){
				for (int j = i; j < basic_cnum - 1; ++j){
					cluster_center[j] = cluster_center[j + 1];
					cluster_node_count[j] = cluster_node_count[j + 1];
				}
				for (int j = 0; j < cluster_index.size(); ++j)
					if (cluster_index[j] > i) cluster_index[j] -= 1;
				basic_cnum--;
			}
		}
	}

	for (int i = 0; i < basic_cnum - 1; ++i)
		for (int j = i + 1; j < basic_cnum; ++j){
			if (sqrt(pow(cluster_center[i][0] - cluster_center[j][0], 2) + pow(cluster_center[i][1] - cluster_center[j][1], 2)) < 1e-5){
				std::cout << "error" << std::endl;
			}
		}

	cluster_center.resize(basic_cnum);
	cluster_node_count.resize(basic_cnum);

	tree_nodes.resize(level_num);
	tree_nodes[0].resize(basic_cnum);
	for (int i = 0; i < basic_cnum; ++i) {
		CLeaf* temp_leaf = new CLeaf();
		temp_leaf->center_pos = cluster_center[i];
		tree_nodes[0][i] = temp_leaf;
		tree_nodes[0][i]->seq_index = i;
	}
	for (int i = 0; i < cluster_index.size(); ++i) {
		CLeaf* leaf_node = dynamic_cast< CLeaf* >(tree_nodes[0][cluster_index[i]]);
		leaf_node->linked_points.push_back(i);
	}
	VtkTriangulation();
}

void LinkedTree::ConstructDirectly(int level_num) {
	ClearData();

	tree_nodes.resize(level_num);
	tree_nodes[0].resize(dataset_->point_pos.size());
	for (int i = 0; i < dataset_->point_pos.size(); ++i) {
		CLeaf* temp_leaf = new CLeaf();
		temp_leaf->center_pos = dataset_->point_pos[i];
		temp_leaf->linked_points.push_back(i);
		tree_nodes[0][i] = temp_leaf;
		tree_nodes[0][i]->seq_index = i;
	}

	if (dataset_->is_structured_data) {
		site_connecting_status.resize(dataset_->point_pos.size());
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			site_connecting_status[i].resize(dataset_->point_pos.size());
			site_connecting_status[i].assign(site_connecting_status[i].size(), false);
		}
		for (int i = 0; i < dataset_->point_pos.size(); ++i) {
			int x = i % dataset_->w;
			int y = i / dataset_->w;
			if (x == dataset_->w - 1 || y == dataset_->h - 1) continue;
			int pre_one = y * dataset_->w + x + 1;
			int pre_two = (y + 1) * dataset_->w + x;
			site_connecting_status[pre_one][i] = true;
			site_connecting_status[i][pre_one] = true;
			site_connecting_status[pre_two][i] = true;
			site_connecting_status[i][pre_two] = true;
		}
	} else {
		VtkTriangulation();
	}
}

void LinkedTree::VtkTriangulation() {
	vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < tree_nodes[0].size(); ++i) {
		CLeaf* leaf_node = dynamic_cast<CLeaf*>(tree_nodes[0][i]);
		points->InsertNextPoint(leaf_node->center_pos[0], leaf_node->center_pos[1], 0);
	}
	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->SetTolerance(1e-30);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

	int site_num = tree_nodes[0].size();
	site_connecting_status.resize(site_num);
	for (int i = 0; i < site_num; ++i) site_connecting_status[i].resize(site_num, false);

	vtkIdTypeArray* idarray = triangle_out->GetPolys()->GetData();
	int pn = triangle_out->GetNumberOfPoints();
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
			CBranch* node_one = dynamic_cast< CBranch* >(tree_nodes[top_level][i]);
			CBranch* node_two = dynamic_cast< CBranch* >(tree_nodes[top_level][j]);
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