#include "gestalt_candidate_set.h"
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

GestaltCandidateSet::GestaltCandidateSet(ScatterPointDataset* data)
	: dataset_(data) {
	this->InitSiteData();
}

GestaltCandidateSet::~GestaltCandidateSet() {

}

void GestaltCandidateSet::InitSiteData() {
	site_num = dataset_->point_pos.size();
	point_site_id.resize(site_num);
	site_point_num.resize(site_num);
	site_center_pos.resize(site_num);
	site_average_value.resize(site_num);
	is_site_labeled.resize(site_num);
	for (int i = 0; i < site_num; ++i) {
		point_site_id[i] = i;
		site_point_num[i] = 1;
		site_center_pos[i].resize(2);
		site_center_pos[i][0] = dataset_->point_pos[i][0];
		site_center_pos[i][1] = dataset_->point_pos[i][1];
		site_average_value[i] = dataset_->point_values[i];
		is_site_labeled[i] = false;
	}

	Triangulation();
}

void GestaltCandidateSet::Triangulation() {
	vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < dataset_->point_pos.size(); ++i) {
		points->InsertNextPoint(dataset_->point_pos[i][0], dataset_->point_pos[i][1], 0);
	}
	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

	int site_num = dataset_->point_pos.size();
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

void GestaltCandidateSet::ExtractGestaltCandidates(float dis_thresh) {
	gestalt_candidates.resize(site_num);

	std::vector< bool > is_reached;
	is_reached.resize(site_num);

	std::vector< float > node_distance;
	node_distance.resize(site_num);

	for (int i = 0; i < site_num; ++i) {
		node_distance.assign(site_num, 1e10);
		is_reached.assign(site_num, false);
		node_distance[i] = 0;

		for (int j = 0; j < site_num; ++j) {
			float min_dist = 1e20;
			int min_index = -1;
			for (int k = 0; k < site_num; ++k)
				if (!is_reached[k] && node_distance[k] < min_dist && node_distance[k] < dis_thresh) {
					min_dist = node_distance[k];
					min_index = k;
				}
			if (min_index == -1) break;
			is_reached[min_index] = true;

			for (int k = 0; k < site_num; ++k)
				if (!is_reached[k] && site_connecting_status[min_index][k]) {
					float temp_dis = sqrt(pow(site_center_pos[min_index][0] - site_center_pos[k][0], 2) + pow(site_center_pos[min_index][1] - site_center_pos[k][1], 2));
					if (node_distance[k] > node_distance[min_index] + temp_dis) node_distance[k] = node_distance[min_index] + temp_dis;
				}
		}
		
		gestalt_candidates[i].clear();
		gestalt_candidates[i].push_back(i);
		for (int j = 0; j < site_num; ++j)
			if (is_reached[j] && j != i) gestalt_candidates[i].push_back(j);
	}
}