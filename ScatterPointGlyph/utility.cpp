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
#include "parallel_coordinate.h"
#include "tour_path_generator.h"

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

void Utility::Sort(std::vector<float>& index_one, std::vector<int>& index_two) {
	for (int i = 0; i < index_one.size() - 1; ++i)
		for (int j = i + 1; j < index_one.size(); ++j)
			if (index_one[i] > index_one[j]){
				float temp_value = index_one[j];
				index_one[j] = index_one[i];
				index_one[i] = temp_value;

				int temp_index = index_two[i];
				index_two[i] = index_two[j];
				index_two[j] = temp_index;
			}
}

void Utility::Triangulation(vector<vector<double>>& pos, vector<vector<bool>>& connecting_status) {
	connecting_status.resize(pos.size());
	for (int i = 0; i < pos.size(); ++i) {
		connecting_status[i].resize(pos.size());
		connecting_status[i].assign(pos.size(), false);
	}
    if (connecting_status.size() <= 2) {
        for (int i = 0; i < connecting_status.size(); ++i)
            for (int j = 0; j < connecting_status[i].size(); ++j)
                connecting_status[i][j] = true;
        return;
    }

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < pos.size(); ++i) {
		points->InsertNextPoint(pos[i][0], pos[i][1], 0);
	}
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->SetTolerance(0.0000001);
	delaunay->SetBoundingTriangulation(false);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

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
	}

    int unconnected_count = 0;
    for (int i = 0; i < connecting_status.size(); ++i){
	    bool is_connecting = false;
	    for (int j = 0; j < connecting_status[i].size(); ++j) is_connecting = is_connecting || connecting_status[i][j];
        if (!is_connecting) {
            int min_index = -1; 
            float min_dis = 1e20;
            for (int j = 0; j < connecting_status[i].size(); ++j) {
                if (j == i) continue;
                float temp_dis;
		        temp_dis = sqrt(pow(pos[i][0] - pos[j][0], 2)
			        + pow(pos[i][1] - pos[j][1], 2));
                if (temp_dis < min_dis) {
                    min_dis = temp_dis;
                    min_index = j;
                }
            }
            if (min_index != -1) {
                connecting_status[i][min_index] = true;
                connecting_status[min_index][i] = true;
                is_connecting = true;
                unconnected_count++;
            }
        }
	    assert(is_connecting);
	}
#ifdef DEBUG_ON
	std::cout << "Unconnected Count: " << unconnected_count << std::endl;
#endif
}

void Utility::GenerateAxisOrder(ParallelDataset* dataset, std::vector<int>& axis_order)
{
	axis_order.resize(dataset->axis_names.size());
	for (int i = 0; i < dataset->axis_names.size(); ++i) axis_order[i] = i;
	std::vector<std::vector<float>> corr_values;
	corr_values.resize(dataset->axis_names.size());
	for (int i = 0; i < dataset->axis_names.size(); ++i) {
		corr_values[i].resize(dataset->axis_names.size(), 0);
	}
	if (dataset->subset_records.size() == 1) {
		for (int i = 0; i < dataset->axis_names.size() - 1; ++i)
			for (int j = i + 1; j < dataset->axis_names.size(); ++j) {
				float temp_corr = 0;
				for (int k = 0; k < dataset->subset_records[0].size(); ++k)
					temp_corr += (dataset->subset_records[0][k]->values[i] - dataset->var_centers[0][i]) * (dataset->subset_records[0][k]->values[j] - dataset->var_centers[0][j]);
				corr_values[i][j] = temp_corr / (dataset->subset_records[0].size() * dataset->var_std_dev[0][i] * dataset->var_std_dev[0][j]);
				corr_values[i][j] = 1.0 - abs(corr_values[i][j]);
				corr_values[j][i] = corr_values[i][j];
			}
		TourPathGenerator::GenerateRoundPath(corr_values, axis_order);
	}
	else {
		std::vector<float> center_values;
		std::vector<float> std_dev_values;
		center_values.resize(dataset->axis_names.size(), 0);
		std_dev_values.resize(dataset->axis_names.size(), 0);
		for (int i = 0; i < dataset->axis_names.size(); ++i) {
			for (int j = 0; j < dataset->subset_records.size(); ++j)
				center_values[i] += dataset->var_centers[j][i];
			center_values[i] /= dataset->subset_records.size();
		}
		for (int i = 0; i < dataset->axis_names.size(); ++i) {
			for (int j = 0; j < dataset->subset_records.size(); ++j)
				std_dev_values[i] += pow(dataset->var_centers[j][i] - center_values[i], 2);
			std_dev_values[i] = sqrt(std_dev_values[i] / dataset->subset_records.size());
		}

		for (int i = 0; i < dataset->axis_names.size() - 1; ++i)
			for (int j = i + 1; j < dataset->axis_names.size(); ++j) {
				float temp_corr = 0;
				for (int k = 0; k < dataset->subset_records.size(); ++k)
					temp_corr += (dataset->var_centers[k][i] - center_values[i]) * (dataset->var_centers[k][j] - center_values[j]);
				corr_values[i][j] = temp_corr / (dataset->subset_records.size() * std_dev_values[i] * std_dev_values[j]);
				corr_values[i][j] = 1.0 - abs(corr_values[i][j]);
				corr_values[j][i] = corr_values[i][j];
			}
		TourPathGenerator::GenerateRoundPath(corr_values, axis_order);
	}
}

void Utility::GridConnection(std::vector<CNode*>& nodes, int w, int h, std::vector<std::vector<bool>>& connecting_status, float& min_edge_length)
{
    connecting_status.resize(nodes.size());
	for (int i = 0; i < nodes.size(); ++i) {
		connecting_status[i].resize(nodes.size());
		connecting_status[i].assign(nodes.size(), false);
	}
    for (int i = 0; i < h; i++) {
        int n1 = i * w;
        int n2 = (i + 1) * w;
        for (int j = 0; j < w; j++) {
            int p1 = n1 + j;
            int p2 = n1 + j + 1;
            int p3 = n2 + j;
            if (i < h - 1) {
                connecting_status[p1][p3] = true;
                connecting_status[p3][p1] = true;
            }
            if (j < w - 1) {
                connecting_status[p1][p2] = true;
                connecting_status[p2][p1] = true;
            }
        }
    }
}

void Utility::GenerateAxisOrder(vector<vector<float>>& values, vector<int>& axis_order) {
    if (values.size() <= 1) return;

    int axis_num = values[0].size();
    int record_num = values.size();

    axis_order.resize(axis_num);
    for (int i = 0; i < axis_num; ++i) axis_order[i] = i;
    if (record_num <= 1) record_num;

    vector<float> average_values;
    vector<float> std_devs;
    average_values.resize(axis_num, 0);
    std_devs.resize(axis_num, 0);
    for (int i = 0; i < record_num; ++i)
        for (int j = 0; j < axis_num; ++j)
            average_values[j] += values[i][j];
    for (int i = 0; i < axis_num; ++i)
        average_values[i] /= record_num;
    for (int i = 0; i < record_num; ++i)
        for (int j = 0; j < axis_num; ++j)
            std_devs[j] += pow(values[i][j] - average_values[j], 2);
    for (int i = 0; i < axis_num; ++i)
        std_devs[i] = sqrt(std_devs[i]);

	std::vector<std::vector<float>> corr_values;
	corr_values.resize(axis_num);
	for (int i = 0; i < axis_num; ++i) {
		corr_values[i].resize(axis_num, 0);
	}
	for (int i = 0; i < axis_num - 1; ++i)
			for (int j = i + 1; j < axis_num; ++j) {
				float temp_corr = 0;
				for (int k = 0; k < record_num; ++k)
					temp_corr += (values[k][i] - average_values[i]) * (values[k][j] - average_values[j]);
				corr_values[i][j] = abs(temp_corr / (std_devs[i] * std_devs[j]));
				corr_values[j][i] = corr_values[i][j];
			}
	TourPathGenerator::GenerateRoundPath(corr_values, axis_order);
}

float Utility::GetAverageDistance(vector<vector<float>>& pos) {
    int edge_count = 0;
    float average_distance = 0;

    if (pos.size() < 5) return GetDirectAverageDistance(pos);

    vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < pos.size(); ++i) {
		points->InsertNextPoint(pos[i][0], pos[i][1], 0);
	}
	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer< vtkDelaunay2D > delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->SetTolerance(0.0000001);
	delaunay->SetBoundingTriangulation(false);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

    if (triangle_out->GetNumberOfPolys() <= 0) return GetDirectAverageDistance(pos);

    vector<float> distances;
	vtkIdTypeArray* idarray = triangle_out->GetPolys()->GetData();
	float min_dis = 1e10;
	for (int i = 0; i < triangle_out->GetNumberOfPolys(); ++i){
		vtkCell* cell = triangle_out->GetCell(i);
		int id1 = cell->GetPointId(0);
		int id2 = cell->GetPointId(1);
		int id3 = cell->GetPointId(2);

		float temp_dis[3];
		temp_dis[0] = sqrt(pow(pos[id1][0] - pos[id2][0], 2)
			+ pow(pos[id1][1] - pos[id2][1], 2));

		temp_dis[1] = sqrt(pow(pos[id1][0] - pos[id3][0], 2)
			+ pow(pos[id1][1] - pos[id3][1], 2));

		temp_dis[2] = sqrt(pow(pos[id3][0] - pos[id2][0], 2)
			+ pow(pos[id3][1] - pos[id2][1], 2));

		average_distance += min(min(temp_dis[0], temp_dis[1]), temp_dis[2]);

        edge_count += 1;
	}
	average_distance /= edge_count;
    return average_distance;
}

bool Utility::CheckInside(vector<float>& path, float x, float y) {
    if (path.size() < 6) return false;

    int count = path.size() / 2;

	if (count < 3) return false;

	const double PI = 3.14159265358979323846;
	double angle = 0;

    double p1[2];
	double p2[2];

	for (int i = 0, k = count - 1; i < count; k = i++) {	
		p1[0] = path[2 * i] - x;
		p1[1] = path[2 * i + 1] - y;

		p2[0] = path[2 * k] - x;
		p2[1] = path[2 * k + 1] - y;

		double radian = atan2(p1[1], p1[0]) - atan2(p2[1], p2[0]);
		radian = abs(radian);
		if (radian > PI)
			radian = 2 * PI - radian;
		angle += radian;
	}
	if (fabs(6.28318530717958647692 - abs(angle)) < 1)
		return true;
	else
		return false;
}

float Utility::GetDirectAverageDistance(vector<vector<float>>& pos) {
    int edge_count = 0;
    float average_distance = 0;

    for (int i = 0; i < pos.size(); ++i) {
        float min_dis = 1e10;
        for (int j = 0; j < pos.size(); ++j) {
            if (j == i) continue;
            float temp_dis = sqrt(pow(pos[i][0] - pos[j][0], 2)
			            + pow(pos[i][1] - pos[j][1], 2));
            min_dis = min(min_dis, temp_dis);
        }
        average_distance += min_dis;
        edge_count++;
    }

    if (edge_count != 0) {
        return average_distance / edge_count;
    }
    else {
        return 1;
    }
}
