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
class QVTKWidget;
class ScatterPointDataset;
class TourPathGenerator;
class vtkPropPicker;

class TransMap : public vtk3DWidget
{
public:
	TransMap(QVTKWidget* parent);

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

	void SetOriginalData(ScatterPointDataset* data);
	void SetData(TransMapData* data);
	void UpdateScale();
	void SetNodeRadius(float r) { this->node_radius = r; }

	void SetBrushSelectionOn();
	void SetBrushSelectionOff();

	int GetSelectedClusterIndex();
	void GetSelectedClusterIndex(std::vector< int >& index);
	void SetSequenceSelectionOn();
	void SetSequenceSelectionOff();

	void SetNodeSelected(int node_id);
	void SetMouseReleased();
	void SetMouseDragmove(int x, int y);

	bool IsMapUpdateNeeded() { return is_map_update_needed;  }

	std::list< CNode* > GetNodeSequence() { return highlight_node_sequence; }

protected:
	TransMap();
	~TransMap();

	std::vector< int > state_vec;
	int state;
	enum WidgetState {
		NORMAL = 0x0,
		HIGHLIGHT_LEVEL_ZERO_NODE,
		HIGHLIGHT_LEVEL_ONE_NODE,
		HIGHLIGHT_TRANSFER_BAND,
		SELECT_PATH_SEQUENCE,
		SELECT_CLUSTER,
		SELECTION_BRUSH_BEGIN,
	};
	
	static void ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata);

	virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();

	void BuildRepresentation();

	std::vector< vtkActor* > level_one_node_glyph_actors;
	std::vector< vtkPolyData* > level_one_node_glyph_polys;

	std::vector< vtkActor* > level_zero_node_glyph_actors;
	std::vector< vtkPolyData* > level_zero_node_glyph_polys;

	std::vector< vtkActor* > trans_glyph_actors;
	std::vector< vtkPolyData* > trans_glyph_polys;
	std::vector< int > trans_edges;

	std::vector< vtkActor* > boundary_glyph_actors;
	std::vector< vtkPolyData* > boundary_glyph_polys;

	QVTKWidget* parent_view;

	vtkActor* highlight_actor;
	vtkPolyDataMapper* hightlight_mapper;
	vtkPolyData* highlight_poly;

	vtkActor* selection_brush_actor;
	vtkPolyDataMapper* selection_brush_mapper;
	vtkPolyData* selection_brush_poly;

	vtkPropPicker* level_one_node_picker;
	vtkPropPicker* level_zero_node_picker;
	vtkPropPicker* trans_picker;
	vtkPropPicker* boundary_picker;

	vtkActor* current_handle;
	CNode* current_node;

	std::list< CNode* > highlight_node_sequence;

	float node_radius;
	bool is_sequence_selection;

	bool is_map_update_needed;

	TourPathGenerator* path_generator_;

private:
	TransMapData* dataset_;
	ScatterPointDataset* scatter_data_;

	void ResizeActors();

	void ConstructActors();
	void HighlightHandle(vtkProp* prop);
	void SelectNode(CNode* node);
	void UpdateHightlightActor();
	void UpdateBrushSelectionResult();
	bool IsInsideSelection(float x, float y);
};

#endif