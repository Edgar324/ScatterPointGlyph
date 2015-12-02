#include "cluster_glyph_layer.h"
#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include "scatter_point_dataset.h"

ClusterGlyphLayer::ClusterGlyphLayer() 
	: BasicGlyphLayer(), dataset_(NULL) {
	radius_range_[0] = 0;
	radius_range_[1] = 0.1;
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

void ClusterGlyphLayer::InitGlyphActor() {
	if (poly_data_ != NULL) {
		poly_data_->Delete();
	}
	poly_data_ = vtkPolyData::New();

	vtkPoints* points = vtkPoints::New();
	poly_data_->SetPoints(points);

	vtkCellArray* poly_array = vtkCellArray::New();
	poly_data_->SetPolys(poly_array);

	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");
	poly_data_->GetPointData()->SetScalars(colors);

	/*vtkCellArray* line_array = vtkCellArray::New();
	poly_data_->SetLines(line_array);
	for (int i = 0; i < line_paras_.size() / 3; ++i) {
	float a = line_paras_[3 * i];
	float b = line_paras_[3 * i + 1];
	float c = line_paras_[3 * i + 2];
	/// ax + by + c = 0
	vtkIdType ids[2];
	ids[0] = poly_data_->GetPoints()->InsertNextPoint(-1 * c / a, 0, 0.005);
	colors->InsertNextTuple3(0, 0, 255);
	ids[1] = poly_data_->GetPoints()->InsertNextPoint(0, -1 * c / b, 0.005);
	colors->InsertNextTuple3(0, 0, 255);
	poly_data_->InsertNextCell(VTK_LINE, 2, ids);
	}*/

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
}

void ClusterGlyphLayer::SetClusterIndex(int cluster_count, std::vector< int >& point_index) {
	if (point_index.size() == 0 || dataset_ == NULL) return;

	InitGlyphActor();

	std::vector< float > x, y;
	x.resize(cluster_count);
	x.assign(cluster_count, 0);
	y.resize(cluster_count);
	y.assign(cluster_count, 0);
	std::vector< int > cluster_point_count;
	cluster_point_count.resize(cluster_count);
	cluster_point_count.assign(cluster_count, 0);
	for (int i = 0; i < point_index.size(); ++i) {
		int cluster_index = point_index[i];
		x[cluster_index] += dataset_->point_pos[i][0];
		y[cluster_index] += dataset_->point_pos[i][1];
		cluster_point_count[cluster_index]++;
	}
	for (int i = 0; i < cluster_count; ++i)
		if (cluster_point_count[i] != 0) {
			x[i] /= cluster_point_count[i];
			y[i] /= cluster_point_count[i];
		} else {
			std::cout << "Empty cluster detected." << std::endl;
			return;
		}

	std::vector< float > cluster_radius;
	cluster_radius.resize(cluster_count, 0);
	for (int i = 0; i < point_index.size(); ++i) {
		int cluster_index = point_index[i];
		cluster_radius[cluster_index] += sqrt(pow(dataset_->point_pos[i][0] - x[cluster_index], 2) + pow(dataset_->point_pos[i][1] - y[cluster_index], 2));
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

	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());

	for (int i = 0; i < cluster_count; ++i) {
		//float r = radius_range_[0] + log((cluster_point_count[i] - min_point_count) / (max_point_count - min_point_count) + 1) / log(2.0) * (radius_range_[1] - radius_range_[0]) * 0.6;
		float r = cluster_radius[i] + radius_range_[0];
		int center_id = poly_data_->GetPoints()->InsertNextPoint(x[i], y[i], 0.01);
		colors->InsertNextTuple3(255, 0, 0);
		int pre_id = poly_data_->GetPoints()->InsertNextPoint(x[i] + r, y[i], 0.01);
		colors->InsertNextTuple3(255, 0, 0);
		vtkIdType cell_ids[3];
		for (int j = 1; j <= 20; ++j) {
			float end_arc = 2 * 3.14159 * j / 20;
			int current_id = poly_data_->GetPoints()->InsertNextPoint(x[i] + r * cos(end_arc), y[i] + r * sin(end_arc), 0.01);
			colors->InsertNextTuple3(255, 0, 0);
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