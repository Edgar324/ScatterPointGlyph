#include "hier_cluster_glyph_layer.h"
#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include "line_extractor.h"
#include "continuity_extractor.h"

HierClusterGlyphLayer::HierClusterGlyphLayer() 
	: BasicGlyphLayer() {

}

HierClusterGlyphLayer::~HierClusterGlyphLayer() {

}

void HierClusterGlyphLayer::SetData(std::vector< float >& pos, std::vector< std::vector< float > >& value) {
	glyph_pos_ = pos;
	glyph_values_ = value;

	highlight_cluster_[0] = -1;
	highlight_cluster_[1] = -1;

	std::vector< float > x;
	std::vector< float > y;
	x.resize(pos.size() / 2);
	y.resize(pos.size() / 2);
	for (int i = 0; i < pos.size() / 2; ++i) {
		x[i] = pos[2 * i];
		y[i] = pos[2 * i + 1];
	}

	
	/*LineExtractor* extractor = new LineExtractor;
	line_paras_.clear();
	extractor->ExtractLines(x, y, line_paras_);*/

	this->UpdateGlyphActor();
}

void HierClusterGlyphLayer::SetHighlighCluster(int cluster_one, int cluster_two) {
	highlight_cluster_[0] = cluster_one;
	highlight_cluster_[1] = cluster_two;

	this->UpdateGlyphActor();
}

void HierClusterGlyphLayer::UpdateGlyphActor() {
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
	/// construct glyph
	for (int i = 0; i < glyph_pos_.size() / 2; ++i) {
		float x = glyph_pos_[i * 2];
		float y = glyph_pos_[i * 2 + 1];
		float r = glyph_values_[i][0] / 2;

		/// insert circle glyph
		int center_id = poly_data_->GetPoints()->InsertNextPoint(x, y, 0.01);
		if (i == highlight_cluster_[0] || i == highlight_cluster_[1])
			colors->InsertNextTuple3(255, 0, 0);
		else
			colors->InsertNextTuple3(0, 255, 0);
		int pre_id = poly_data_->GetPoints()->InsertNextPoint(x + r, y, 0.01);
		if (i == highlight_cluster_[0] || i == highlight_cluster_[1])
			colors->InsertNextTuple3(255, 0, 0);
		else
			colors->InsertNextTuple3(0, 255, 0);
		vtkIdType cell_ids[3];
		for (int j = 1; j <= 20; ++j) {
			float end_arc = 2 * 3.14159 * j / 20;
			int current_id = poly_data_->GetPoints()->InsertNextPoint(x + r * cos(end_arc), y + r * sin(end_arc), 0.01);
			if (i == highlight_cluster_[0] || i == highlight_cluster_[1])
				colors->InsertNextTuple3(255, 0, 0);
			else
				colors->InsertNextTuple3(0, 255, 0);
			cell_ids[0] = center_id;
			cell_ids[1] = pre_id;
			cell_ids[2] = current_id;
			poly_data_->InsertNextCell(VTK_TRIANGLE, 3, cell_ids);

			pre_id = current_id;
		}
	}

	vtkCellArray* line_array = vtkCellArray::New();
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
	}

	poly_data_->GetPointData()->SetScalars(colors);

	mapper_->SetInputData(poly_data_);
	actor_->SetMapper(mapper_);
}