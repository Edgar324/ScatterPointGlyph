#ifndef FIELD_RENDERING_LAYER_H_
#define FIELD_RENDERING_LAYER_H_

#include <vtk3DWidget.h>
#include <vtkScalarBarActor.h>
#include <vector>

class vtkActor;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkLookupTable;
class ScatterPointDataset;

class FieldRenderingLayer : public vtk3DWidget
{
public:
	FieldRenderingLayer();
	~FieldRenderingLayer();

	vtkTypeMacro(FieldRenderingLayer, vtk3DWidget);
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

	void SetFieldData(ScatterPointDataset* data, int var_index);
	void SetColorScalar(std::vector< float >& values, std::vector< double >& rgb);

private:
	vtkPolyData* field_polydata_;
	vtkPolyDataMapper* field_poly_mapper_;
	vtkActor* field_actor_;

	vtkPolyData* point_polydata_;
	vtkPolyDataMapper* point_mapper_;
	vtkActor* point_actor_;

	vtkScalarBarActor* color_bar_actor_;
	vtkLookupTable* lookup_table_;

	ScatterPointDataset* dataset_;
	int var_index_;
	bool is_segment_probability_;

	void ConstructFieldActors();
};

#endif