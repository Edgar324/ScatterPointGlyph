#ifndef DISTANCE_MATRIX_LAYER_H_
#define DISTANCE_MATRIX_LAYER_H_

#include <vtk3DWidget.h>
#include <vector>
#include <string>

class vtkActor;
class vtkTextWidget;
class vtkTextActor;
class vtkPolyData;
class vtkPolyDataMapper;

class DistanceMatrixLayer : public vtk3DWidget
{
public:
	DistanceMatrixLayer();
	~DistanceMatrixLayer();

	vtkTypeMacro(DistanceMatrixLayer, vtk3DWidget);
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

	void SetData(std::vector< std::string >& names, std::vector< std::vector< float > >& distance);
	void BuildRepresentation();

private:
	vtkActor* actor_;
	vtkPolyDataMapper* mapper_;
	vtkPolyData* poly_data_;

	vtkTextActor* title_actor_;

	std::vector< vtkTextActor* > hor_label_actors_, ver_label_actors_;

	std::vector< std::string > label_names_;
	std::vector< std::vector< float > > label_distance_;
};

#endif