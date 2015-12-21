#ifndef TRANSMAP_H_
#define TRANSMAP_H_

#include <vtk3DWidget.h>
#include <vtkLineWidget.h>
#include <vector>

class vtkActor;
class vtkPolyDataMapper;
class vtkCellPicker;
class vtkPolyData;
class vtkProperty;
class TransMapData;

class TransMap : public vtk3DWidget
{
public:
	static TransMap* New();

	vtkTypeMacro(TransMap, vtk3DWidget);
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

	void SetData(TransMapData* data);
	void SetNodeRadius(float r) { this->node_radius = r; }

protected:
	TransMap();
	~TransMap();

	int state;
	enum WidgetState {
		START = 0x0,
	};
	
	static void ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata);

	virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();

	void BuildRepresentation();

	std::vector< vtkActor* > node_glyph_actors;
	std::vector< vtkPolyDataMapper* > node_glyph_mappers;
	std::vector< vtkPolyData* > node_glyph_polys;

	std::vector< vtkActor* > trans_glyph_actors;
	std::vector< vtkPolyDataMapper* > trans_glyph_mappers;
	std::vector< vtkPolyData* > trans_glyph_polys;

	std::vector< vtkActor* > boundary_glyph_actors;
	std::vector< vtkPolyDataMapper* > boundary_glyph_mappers;
	std::vector< vtkPolyData* > boundary_glyph_polys;

	vtkCellPicker* node_picker;
	vtkCellPicker* trans_picker;
	vtkCellPicker* boundary_picker;
	vtkActor* current_handle;

	float node_radius;

private:
	TransMapData* dataset_;

	void ClearActors();
	void InitActors();
};

#endif