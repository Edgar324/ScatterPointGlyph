#ifndef POINT_RENDERING_LAYER_H_
#define POINT_RENDERING_LAYER_H_

#include <vector>
#include <vtk3DWidget.h>
#include <QtGui/QColor>

class ScatterPointDataset;
class vtkActor;
class vtkPolyData;
class vtkPolyDataMapper;

class PointRenderingLayer : public vtk3DWidget
{
public:
	PointRenderingLayer();
	~PointRenderingLayer();

	vtkTypeMacro(PointRenderingLayer, vtk3DWidget);
	void PrintSelf(ostream& os, vtkIndent indent) {}

	virtual void SetEnabled(int);
	virtual void PlaceWidget(double bounds[6]) {}
	void PlaceWidget() {
		this->Superclass::PlaceWidget();
	}
	void PlaceWidget(double xmin, double xmax, double ymin, double ymax,
		double zmin, double zmax) {
		this->Superclass::PlaceWidget(xmin, xmax, ymin, ymax, zmin, zmax);
	}

	void SetData(ScatterPointDataset* data);
	void SetPointValue(std::vector< float >& values);

	void SetClusterIndex(int cluster_count, std::vector< int >& point_index, std::vector< QColor >& colors);
	void SetHighlightCluster(int index);
	void SetHighlightClusters(std::vector< int >& index);

	void SetCategoryOn();
	void SetCategoryOff();

	std::vector< int >& GetClusterColor() { return cluster_color_; }

private:
	int cluster_count_;
	std::vector< int > point_index_;
	std::vector< float > point_values_;
	std::vector< int > cluster_color_;

	bool is_category_on_;
	std::vector< int > current_selection_index_;

	ScatterPointDataset* dataset_;

	vtkActor* actor_;
	vtkPolyDataMapper* mapper_;
	vtkPolyData* poly_data_;

	void UpdatePointGlyph();
};

#endif