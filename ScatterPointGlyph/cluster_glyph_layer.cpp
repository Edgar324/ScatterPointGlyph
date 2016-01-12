#include "cluster_glyph_layer.h"
#include <algorithm>
#include <vector> 
#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include "scatter_point_dataset.h"

ClusterGlyphLayer::ClusterGlyphLayer() 
	: BasicGlyphLayer(), dataset_(NULL) {
	radius_range_[0] = 0;
	radius_range_[1] = 0.1;
	current_cluster_count_ = -1;
}

ClusterGlyphLayer::~ClusterGlyphLayer() {

}

void ClusterGlyphLayer::SetData(ScatterPointDataset* data) {
	dataset_ = data;

	this->InitGlyphActor();
}

void ClusterGlyphLayer::SetRadiusRange(float maxR, float minR){
	radius_range_[0] = minR;
	radius_range_[1] = maxR;
}

void ClusterGlyphLayer::Select(float x, float y) {
	int nearest_index = -1;
	float min_dis = 1e10;
	for (int i = 0; i < center_x_.size(); ++i) {
		if (abs(x - center_x_[i]) + abs(y - center_y_[i]) < min_dis) {
			min_dis = abs(x - center_x_[i]) + abs(y - center_y_[i]);
			nearest_index = i;
		}
	}
	if (nearest_index != -1) {
		is_cluster_selected_[nearest_index] = !is_cluster_selected_[nearest_index];
	}
}

void ClusterGlyphLayer::InitGlyphActor() {
	if (poly_data_ == NULL) {
		poly_data_ = vtkPolyData::New();
	} else {
		poly_data_->Initialize();
	}
	

	vtkPoints* points = vtkPoints::New();
	poly_data_->SetPoints(points);

	vtkCellArray* poly_array = vtkCellArray::New();
	poly_data_->SetPolys(poly_array);

	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");
	poly_data_->GetPointData()->SetScalars(colors);

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
}

void ClusterGlyphLayer::SetClusterIndex(int cluster_count, std::vector< int >& point_index) {
	if (point_index.size() == 0 || dataset_ == NULL || current_cluster_count_ == cluster_count) return;
	current_cluster_count_ = cluster_count;

	InitGlyphActor();

	is_cluster_selected_.resize(cluster_count);
	is_cluster_selected_.assign(cluster_count, false);

	center_x_.resize(cluster_count);
	center_x_.assign(cluster_count, 0);
	center_y_.resize(cluster_count);
	center_y_.assign(cluster_count, 0);
	std::vector< std::vector< float > > average_values;
	average_values.resize(cluster_count);
	for (int i = 0; i < cluster_count; ++i) average_values[i].resize(dataset_->var_weights.size(), 0);

	std::vector< int > cluster_point_count;
	cluster_point_count.resize(cluster_count);
	cluster_point_count.assign(cluster_count, 0);
	for (int i = 0; i < point_index.size(); ++i) {
		int cluster_index = point_index[i];
		if (cluster_index < 0) continue;
		center_x_[cluster_index] += dataset_->original_point_pos[i][0];
		center_y_[cluster_index] += dataset_->original_point_pos[i][1];
		for (int j = 0; j < dataset_->var_weights.size(); ++j)
			average_values[cluster_index][j] += (dataset_->original_point_values[i][j] - dataset_->original_value_ranges[j][0]) / (dataset_->original_value_ranges[j][1] - dataset_->original_value_ranges[j][0]);
		cluster_point_count[cluster_index]++;
	}
	for (int i = 0; i < cluster_count; ++i)
		if (cluster_point_count[i] != 0) {
			center_x_[i] /= cluster_point_count[i];
			center_y_[i] /= cluster_point_count[i];
			for (int j = 0; j < dataset_->var_weights.size(); ++j)
				average_values[i][j] /= cluster_point_count[i];
		} else {
			std::cout << "Empty cluster detected." << std::endl;
			return;
		}

	std::vector< int > new_count;
	new_count = cluster_point_count;
	std::sort(new_count.data(), new_count.data() + new_count.size());
	int thresh = -1;
	for (int i = 0; i < cluster_point_count.size(); ++i)
		if (cluster_point_count[i] > thresh) thresh = cluster_point_count[i];
	thresh = (int)(thresh * 0.1);

	std::vector< float > cluster_radius;
	cluster_radius.resize(cluster_count, 0);
	for (int i = 0; i < point_index.size(); ++i) {
		int cluster_index = point_index[i];
		if (cluster_index < 0) continue;
		cluster_radius[cluster_index] += sqrt(pow(dataset_->original_point_pos[i][0] - center_x_[cluster_index], 2) + pow(dataset_->original_point_pos[i][1] - center_y_[cluster_index], 2));
	}
	for (int i = 0; i < cluster_count; ++i) {
		cluster_radius[i] /= cluster_point_count[i];
	}

	float max_point_count = 0;
	float min_point_count = 10000000;
	for (int i = 0; i < cluster_count; ++i) {
		if (cluster_point_count[i] > max_point_count) max_point_count = cluster_point_count[i];
		if (cluster_point_count[i] < min_point_count) min_point_count = cluster_point_count[i];
	}

	int red[] = { 166, 223, 128 };
	int green[] = { 97, 194, 205 };
	int blue[] = { 26, 125, 193 };

	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());

	for (int i = 0; i < cluster_count; ++i) {
		//float r = radius_range_[0] + log((cluster_point_count[i] - min_point_count) / (max_point_count - min_point_count) + 1) / log(2.0) * (radius_range_[1] - radius_range_[0]) * 0.6;
		float r = radius_range_[0] * 0.8;
		float base_r = radius_range_[0] * 0.2;
		if (cluster_point_count[i] <= thresh) r = 0.5 * radius_range_[0];
		int center_id = poly_data_->GetPoints()->InsertNextPoint(center_x_[i], center_y_[i], 0.000001);
		/*if (cluster_point_count[i] <= thresh)
			colors->InsertNextTuple3(0, 255, 0);
		else*/
		colors->InsertNextTuple3(red[0], green[0], blue[0]);
		int pre_id = poly_data_->GetPoints()->InsertNextPoint(center_x_[i] + r * average_values[i][0] + base_r, center_y_[i], 0.000001);
		/*if (cluster_point_count[i] <= thresh)
			colors->InsertNextTuple3(0, 255, 0);
		else*/
		colors->InsertNextTuple3(red[0], green[0], blue[0]);
		vtkIdType cell_ids[3];
		for (int j = 1; j <= 20 * dataset_->var_weights.size(); ++j) {
			int var_index = (j - 1) / 20;
			float end_arc = 2 * 3.14159 * j / (20 * dataset_->var_weights.size());
			int current_id = poly_data_->GetPoints()->InsertNextPoint(center_x_[i] + r * cos(end_arc) * average_values[i][var_index] + base_r * cos(end_arc), center_y_[i] + r * sin(end_arc) * average_values[i][var_index] + base_r * sin(end_arc), 0.000001);
			/*if (cluster_point_count[i] <= thresh)
				colors->InsertNextTuple3(0, 255, 0);
			else*/
			colors->InsertNextTuple3(red[var_index], green[var_index], blue[var_index]);
			cell_ids[0] = center_id;
			cell_ids[1] = pre_id;
			cell_ids[2] = current_id;
			poly_data_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
			pre_id = current_id;
		}
	}
	colors->Modified();
	actor_->Modified();
}

void ClusterGlyphLayer::ClearGlyph() {

}