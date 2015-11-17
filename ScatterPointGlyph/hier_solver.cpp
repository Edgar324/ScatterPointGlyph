#include "hier_solver.h"
#include <time.h>
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

HierSolver::HierSolver() 
	: current_cluster_num_(-1), expected_cluster_num_(-1), dis_scale_(0.3) {

}

HierSolver::~HierSolver() {

}

void HierSolver::SetData(ScatterPointDataset* data) {
	dataset_ = data;

	final_label_.resize(dataset_->point_pos.size());
	cluster_center_.resize(dataset_->point_pos.size());
	for (size_t i = 0; i < final_label_.size(); ++i) {
		final_label_[i] = i;
		cluster_center_[i].resize(3);
		cluster_center_[i][0] = dataset_->point_pos[i][0];
		cluster_center_[i][1] = dataset_->point_pos[i][1];
		cluster_center_[i][2] = dataset_->point_values[i];
	}

	cluster_node_count_.resize(dataset_->point_pos.size());
	cluster_node_count_.assign(cluster_node_count_.size(), 1);
	is_label_used_.resize(final_label_.size());
	is_label_used_.assign(is_label_used_.size(), true);
	current_cluster_num_ = final_label_.size();

	Triangulation();
}

void HierSolver::SetExpectedClusterNum(int num) {
	expected_cluster_num_ = num;
}

void HierSolver::Triangulation() {
	vtkSmartPointer< vtkPoints > points = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i < dataset_->point_pos.size(); ++i) {
		points->InsertNextPoint(dataset_->point_pos[i][0], dataset_->point_pos[i][1], 0);
	}
	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

	int cluster_num = dataset_->point_pos.size();
	cluster_connecting_status_.resize(cluster_num);
	for (int i = 0; i < cluster_num; ++i) cluster_connecting_status_[i].resize(cluster_num, false);

	vtkIdTypeArray* idarray = triangle_out->GetPolys()->GetData();
	for (int i = 0; i < triangle_out->GetNumberOfPolys(); ++i){
		double* temp = idarray->GetTuple(i * 3);
		int id1 = (int)(temp[0]);
		temp = idarray->GetTuple(i * 3 + 1);
		int id2 = (int)(temp[0]);
		temp = idarray->GetTuple(i * 3 + 2);
		int id3 = (int)(temp[0]);
		cluster_connecting_status_[id1][id2] = true;
		cluster_connecting_status_[id2][id1] = true;
		cluster_connecting_status_[id2][id3] = true;
		cluster_connecting_status_[id3][id2] = true;
		cluster_connecting_status_[id1][id3] = true;
		cluster_connecting_status_[id3][id1] = true;
	}

	for (int i = 0; i < cluster_connecting_status_.size(); ++i){
		bool is_connecting = false;
		for (int j = 0; j < cluster_connecting_status_[i].size(); ++j) is_connecting = is_connecting || cluster_connecting_status_[i][j];
		assert(is_connecting);
	}
}

void HierSolver::run() {
	while (current_cluster_num_ > expected_cluster_num_) {
		int temp_index[2];
		float min_dis = 1e10;
		for (int i = 0; i < cluster_center_.size() - 1; ++i)
			if (is_label_used_[i]) {
				for (int j = i + 1; j < cluster_center_.size(); ++j)
					if (is_label_used_[j] && cluster_connecting_status_[i][j]) {
						float temp_dis = sqrt(pow(cluster_center_[i][0] - cluster_center_[j][0], 2) + pow(cluster_center_[i][1] - cluster_center_[j][1], 2)) * dis_scale_;
						temp_dis += abs(cluster_center_[i][2] - cluster_center_[j][2]) * (1.0 - dis_scale_);
						if (temp_dis < min_dis) {
							min_dis = temp_dis;
							temp_index[0] = i;
							temp_index[1] = j;
						}
					}
			}

		// update center and cluster index
		for (int i = 0; i < 3; ++i)
			cluster_center_[temp_index[0]][i] = cluster_center_[temp_index[0]][i] * cluster_node_count_[temp_index[0]]
				+ cluster_center_[temp_index[1]][i] * cluster_node_count_[temp_index[1]];

		cluster_node_count_[temp_index[0]] += cluster_node_count_[temp_index[1]];
		for (int i = 0; i < final_label_.size(); ++i)
			if (final_label_[i] == temp_index[1]) final_label_[i] = temp_index[0];

		for (int i = 0; i < 3; ++i) 
			cluster_center_[temp_index[0]][i] /= cluster_node_count_[temp_index[0]];

		is_label_used_[temp_index[1]] = false;
		for (int i = 0; i < cluster_connecting_status_.size(); ++i)
			cluster_connecting_status_[temp_index[0]][i] = cluster_connecting_status_[temp_index[0]][i] || cluster_connecting_status_[temp_index[1]][i];

		current_cluster_num_--;
	}

	// update final label
	final_label_count_ = expected_cluster_num_;
	int accu_count = 0;
	std::vector< int > assigned_labels;
	assigned_labels.resize(is_label_used_.size(), -1);
	for (int i = 0; i < is_label_used_.size(); ++i) {
		if (is_label_used_[i]) {
			assigned_labels[i] = accu_count;
			accu_count++;
		}
	}
	final_label_.resize(dataset_->point_pos.size());
	for (int i = 0; i < final_label_.size(); ++i) {
		final_label_[i] = assigned_labels[final_label_[i]];
	}
}