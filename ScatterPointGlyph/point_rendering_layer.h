#ifndef POINT_RENDERING_LAYER_H_
#define POINT_RENDERING_LAYER_H_

#include <vector>
#include <vtk3DWidget.h>
#include <QtGui/QColor>

class ScatterPointDataset;
class vtkActor;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkScalarBarActor;
class vtkLookupTable;
class vtkRenderer;

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

    void SetMapEnabled(bool enabled);
    void SetScatterPointEnabled(bool enabled);

    void SetColorBarRenderer(vtkRenderer* renderer);

	std::vector< int >& GetClusterColor() { return cluster_color_; }

	void ClearView();

private:
	int cluster_count_;
	std::vector< int > point_index_;
	std::vector< float > point_values_;
	std::vector< int > cluster_color_;

	bool is_category_on_;
	std::vector< int > current_selection_index_;

	ScatterPointDataset* dataset_;

    vtkRenderer* color_bar_renderer_;

	vtkActor* actor_;
	vtkPolyDataMapper* mapper_;
	vtkPolyData* poly_data_;

    vtkActor* map_actor_;
	vtkPolyDataMapper* map_mapper_;
	vtkPolyData* map_poly_data_;

    vtkScalarBarActor* bar_actor_;
    vtkLookupTable* scalar_lookup_table_;

	void UpdatePointGlyph();
    void UpdateValueMapping();
    void LoadMap(const char* file_name, float start_x, float end_x, float start_y, float end_y);
};

#endif