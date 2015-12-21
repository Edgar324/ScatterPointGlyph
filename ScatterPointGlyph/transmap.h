#ifndef TRANSMAP_H_
#define TRANSMAP_H_

#include <vtk3DWidget.h>
#include <QtCore/QObject>
#include <vector>

class vtkActor;
class vtkPolyDataMapper;
class vtkCellPicker;
class vtkPolyData;
class vtkProperty;
class TransMapData;
class vtkProp;
class CNode;

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

	int GetSelectedClusterIndex();

protected:
	TransMap();
	~TransMap();

	int state;
	enum WidgetState {
		NORMAL = 0x0,
		HIGHLIGHT_LEVEL_ZERO_NODE,
		HIGHLIGHT_LEVEL_ONE_NODE,
		HIGHLIGHT_TRANSFER_BAND,
		SELECT_PATH_SEQUENCE,
	};
	
	static void ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata);

	virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();

	void BuildRepresentation();

	std::vector< vtkActor* > level_one_node_glyph_actors;
	std::vector< vtkPolyDataMapper* > node_glyph_mappers;
	std::vector< vtkPolyData* > node_glyph_polys;

	std::vector< vtkActor* > level_zero_node_glyph_actors;

	std::vector< vtkActor* > trans_glyph_actors;
	std::vector< vtkPolyDataMapper* > trans_glyph_mappers;
	std::vector< vtkPolyData* > trans_glyph_polys;

	std::vector< vtkActor* > boundary_glyph_actors;
	std::vector< vtkPolyDataMapper* > boundary_glyph_mappers;
	std::vector< vtkPolyData* > boundary_glyph_polys;

	vtkActor* highlight_actor;
	vtkPolyDataMapper* hightlight_mapper;
	vtkPolyData* highlight_poly;

	vtkCellPicker* level_one_node_picker;
	vtkCellPicker* level_zero_node_picker;
	vtkCellPicker* trans_picker;
	vtkCellPicker* boundary_picker;

	vtkActor* current_handle;
	CNode* current_node;

	float node_radius;

private:
	TransMapData* dataset_;

	void ClearActors();

	void ConstructActors();
	void HighlightHandle(vtkProp* prop);
};

#endif