#include "cluster_glyph_layer.h"
#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include "line_extractor.h"
#include "continuity_extractor.h"

ClusterGlyphLayer::ClusterGlyphLayer() 
	: BasicGlyphLayer() {

}

ClusterGlyphLayer::~ClusterGlyphLayer() {

}

void ClusterGlyphLayer::SetData(std::vector< std::vector< float > >& pos, std::vector< std::vector< float > >& value) {
	glyph_pos_ = pos;
	glyph_values_ = value;

	highlight_cluster_[0] = -1;
	highlight_cluster_[1] = -1;

	
	/*LineExtractor* extractor = new LineExtractor;
	line_paras_.clear();
	extractor->ExtractLines(x, y, line_paras_);*/

	this->UpdateGlyphActor();
}

void ClusterGlyphLayer::SetHighlighCluster(int cluster_one, int cluster_two) {
	highlight_cluster_[0] = cluster_one;
	highlight_cluster_[1] = cluster_two;

	this->UpdateGlyphActor();
}

void ClusterGlyphLayer::AddGestaltGlyph(std::vector< int >& point_index) {
	if (point_index.size() == 0) return;

	/*if (poly_data_ != NULL) {
		poly_data_->Delete();
		poly_data_ = vtkPolyData::New();
	}

	vtkPoints* points = vtkPoints::New();
	poly_data_->SetPoints(points);
	vtkCellArray* poly_array = vtkCellArray::New();
	poly_data_->SetPolys(poly_array);

	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");*/
	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::SafeDownCast(poly_data_->GetPointData()->GetScalars());

	float x = 0;
	float y = 0;
	for (int i = 0; i < point_index.size(); ++i) {
		x += glyph_pos_[point_index[i]][0];
		y += glyph_pos_[point_index[i]][1];
	}
	x /= point_index.size();
	y /= point_index.size();

	/// construct glyph
	float r = 0.005 * log(point_index.size() + 1) + 0.001;

	/// insert circle glyph
	int center_id = poly_data_->GetPoints()->InsertNextPoint(x, y, 0.01);
	colors->InsertNextTuple3(255, 0, 0);
	int pre_id = poly_data_->GetPoints()->InsertNextPoint(x + r, y, 0.01);
	colors->InsertNextTuple3(255, 0, 0);
	vtkIdType cell_ids[3];
	for (int j = 1; j <= 20; ++j) {
		float end_arc = 2 * 3.14159 * j / 20;
		int current_id = poly_data_->GetPoints()->InsertNextPoint(x + r * cos(end_arc), y + r * sin(end_arc), 0.01);
		colors->InsertNextTuple3(255, 0, 0);
		cell_ids[0] = center_id;
		cell_ids[1] = pre_id;
		cell_ids[2] = current_id;
		poly_data_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);
		pre_id = current_id;
	}

	poly_data_->GetPointData()->SetScalars(colors);

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
	actor_->Modified();
}

void ClusterGlyphLayer::ClearGlyph() {

}

void ClusterGlyphLayer::UpdateGlyphActor() {
	if (poly_data_ != NULL) {
		poly_data_->Delete();
	}
	poly_data_ = vtkPolyData::New();

	vtkPoints* points = vtkPoints::New();
	poly_data_->SetPoints(points);

	vtkCellArray* poly_array = vtkCellArray::New();
	poly_data_->SetPolys(poly_array);

	vtkCellArray* line_array = vtkCellArray::New();
	poly_data_->SetLines(line_array);

	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");
	poly_data_->GetPointData()->SetScalars(colors);
	
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
	}

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
}