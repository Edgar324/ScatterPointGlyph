#ifndef TRANSMAP_H_
#define TRANSMAP_H_

#include <vtk3DWidget.h>
#include <QtCore/QObject>
#include <QtCore/QTimer>
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
	void SetAxisOrder(std::vector<int>& order);
    void SetIndicatorRenderer(vtkRenderer* renderer);

	void ClearView();

	std::list<CNode*> GetNodeSequence() { return highlight_node_sequence; }

	void HighlightVar(int var_index);
	void ShowMinimumSpanningTree(bool enabled);
	void ShowVarTrend(int var_index);

	void OnNodeSelected(int node_id);
	void OnMouseReleased(bool is_left_button = true);
	void OnMouseMove(int x, int y);

    void GetFocusVarIndex(std::vector<int>& focus_index);
	int GetSelectedClusterIndex();
	void GetSelectedClusterIndex(std::vector<int>& index);
	void GetSelectedClusterIds(std::vector<int>& ids);
	void GetSelectedClusterNodes(std::vector<CNode*>& nodes);

    void SetDensityMapVisibility(bool visibility);
    void UpdateDensityActor(std::vector<QColor >& node_colors);

	void ForceFocusCenter();
	void MoveViewToFocus(float scale);

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

	std::vector<vtkActor* > node_glyph_actors;
	std::vector<vtkPolyData* > node_glyph_polys;

    std::vector<int> linked_edges;
    vtkActor* linked_glyph_actors;
	vtkPolyDataMapper* linked_glyph_mapper;
	vtkPolyData* linked_glyph_polys;

    std::vector<int> trans_edges;
	vtkActor* trans_glyph_actors;
	vtkPolyDataMapper* trans_glyph_mapper;
	vtkPolyData* trans_glyph_polys;

	vtkActor* highlight_actor;
	vtkPolyDataMapper* hightlight_mapper;
	vtkPolyData* highlight_poly;
	std::vector<vtkTextActor3D* > seqence_text_actors;

	vtkActor* selection_brush_actor;
	vtkPolyDataMapper* selection_brush_mapper;
	vtkPolyData* selection_brush_poly;

    vtkActor* density_actor_;
	vtkPolyDataMapper* density_mapper_;
	vtkPolyData* density_poly_;

    vtkActor* indicator_actor_;
	vtkPolyDataMapper* indicator_mapper_;
	vtkPolyData* indicator_poly_;
    vtkTextActor3D* indicator_text_;

    vtkActor* var_icon_actor_;
    vtkPolyDataMapper* var_icon_mapper_;
    vtkPolyData* var_icon_poly_;
    std::vector<vtkTextActor3D* > variable_color_text_;

	vtkPropPicker* node_picker;
	vtkTooltipItem* tool_tip_item_;

	std::list<CNode*> highlight_node_sequence;

	TourPathGenerator* path_generator_;

private:
	TransMapData* dataset_;
	ScatterPointDataset* scatter_data_;
	float node_radius_;
	std::vector<int> axis_order_;

    vtkRenderer* glyph_indicator_renderer_;

	CNode* current_node_;

	bool is_mst_fixed_ = false;
	bool is_var_trend_fixed_ = false;
	int var_trend_index_ = -1;
    bool is_hovering_ = false;
    bool is_densitymap_shown_ = true;

    int current_highlight_var_index_;

	std::vector<int> focus_var_index_;
	bool is_focus_var_fixed_ = false;

	double focus_pos_[2];
	double current_pos_[3];
	double eye_pos_[3];

	void UpdateNodeActors();
	void GenerateTransEdgeFromHighlight();
	void UpdateTransEdgeActor();
	void UpdateHightlightActor();
	void UpdateIconActor();

	void OnNodeSelected(CNode* node);

	CNode* GetSelectedNode(int x, int y);
	int GetLevelOneNodeIndex(CNode* node);
	int GetClusterNodeIndex(CNode* node);
	bool IsLevelOneNode(CNode* node);

	bool IsInsideSelection(float x, float y);
	void Sort(std::vector<int>& index_one, std::vector<int>& index_two);
};

#endif