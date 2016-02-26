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
class vtkTextActor3D;
class vtkTooltipItem;

class ScatterPointView;

class TransMap : public vtk3DWidget
{
public:
	TransMap(ScatterPointView* parent);

	vtkTypeMacro(TransMap, vtk3DWidget);
	void PrintSelf(ostream& os, vtkIndent indent) {}

	enum WidgetState {
		NORMAL = 0x0,
		SELECT_SINGLE_CLUSTER,
		SELECT_MULTI_CLUSTERS,
		SELECT_MINIMUM_PATH,
		SELECT_BRUSHED_PATH_SEQUENCE,
		SELECT_BRUSHED_CLUSTERS,
	};

	virtual void SetEnabled(int);
	virtual void PlaceWidget(double bounds[6]) {}
	void PlaceWidget() {
		this->Superclass::PlaceWidget();
	}
	void PlaceWidget(double xmin, double xmax, double ymin, double ymax,
		double zmin, double zmax) {
		this->Superclass::PlaceWidget(xmin, xmax, ymin, ymax, zmin, zmax);
	}

	void SetData(ScatterPointDataset* ori_data, TransMapData* data);
	void SetNodeRadius(float r);
	void SetInteractionState(WidgetState s);
	void SetAxisOrder(std::vector< int >& order);

	std::list< CNode* > GetNodeSequence() { return highlight_node_sequence; }

	void HighlightVar(int var_index);
	void ShowMinimumSpanningTree(bool enabled);
	void ShowVarTrend(int var_index);

	void OnNodeSelected(int node_id);
	void OnMouseReleased(bool is_left_button = true);
	void OnMouseMove(int x, int y);

	int GetSelectedClusterIndex();
	void GetSelectedClusterIndex(std::vector< int >& index);
	void GetSelectedClusterIds(std::vector< int >& ids);
	void GetSelectedClusterNodes(std::vector< CNode* >& nodes);

protected:
	TransMap();
	~TransMap();

	int state;
	
	static void ProcessEvents(vtkObject* object, unsigned long event, void* clientdata, void* calldata);

	virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp() {}

	void BuildRepresentation();

	ScatterPointView* parent_view;

	std::vector< vtkActor* > node_glyph_actors;
	std::vector< vtkPolyData* > node_glyph_polys;
	std::vector< vtkTextActor3D* > point_num_text_actors;

	vtkActor* trans_glyph_actors;
	vtkPolyDataMapper* trans_glyph_mapper;
	vtkPolyData* trans_glyph_polys;

	vtkActor* highlight_actor;
	vtkPolyDataMapper* hightlight_mapper;
	vtkPolyData* highlight_poly;
	std::vector< vtkTextActor3D* > seqence_text_actors;

	vtkActor* selection_brush_actor;
	vtkPolyDataMapper* selection_brush_mapper;
	vtkPolyData* selection_brush_poly;

	vtkPropPicker* node_picker;
	vtkTooltipItem* tool_tip_item_;

	std::vector< int > trans_edges;
	std::list< CNode* > highlight_node_sequence;

	TourPathGenerator* path_generator_;

private:
	TransMapData* dataset_;
	ScatterPointDataset* scatter_data_;
	float node_radius_;
	std::vector< int > axis_order_;

	CNode* current_node_;

	bool is_mst_fixed_ = false;
	bool is_var_trend_fixed_ = false;
	int var_trend_index_ = -1;
	int highlight_var_index_ = -1;
	bool is_highlight_var_fixed_ = false;

	void UpdateNodeActors();
	void GenerateTransEdgeFromHighlight();
	void UpdateTransEdgeActor();
	void UpdateHightlightActor();
	void UpdateTipActor();

	void OnNodeSelected(CNode* node);

	CNode* GetSelectedNode(int x, int y);
	int GetLevelOneNodeIndex(CNode* node);
	int GetClusterNodeIndex(CNode* node);
	bool IsLevelOneNode(CNode* node);

	bool IsInsideSelection(float x, float y);
	void Sort(std::vector< int >& index_one, std::vector< int >& index_two);
};

#endif