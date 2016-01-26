#ifndef TREE_MAP_VIEW_H_
#define TREE_MAP_VIEW_H_

#include <QtWidgets/QGraphicsView>
#include "tree_common.h"

class QGraphicsTextItem;
class TreeMapItem;
class VariableItem;

class TreeMapView : public QGraphicsView
{
	Q_OBJECT

public:
	TreeMapView();
	~TreeMapView();

	void SetData(CNode* root, int var_num, std::vector< CNode* >& selected_nodes, int selected_count, std::vector< int >& order, std::vector< QString >& names);
	void SetHighlightVarIndex(int index);
	void SetTreeMapVisible(bool visible);
	void SetTableLensVisible(bool visible);
	bool IsVisible() { return is_table_lens_visible_ || is_treemap_visible_; }

signals:
	void NodeSelected(int node_id);
	void HighlightVarChanged(int var_index);

private:
	CNode* root_node_;
	int var_num_;
	std::vector< int > var_order_;

	bool is_treemap_visible_;
	bool is_table_lens_visible_;
	int highlight_var_index_;

	QGraphicsScene* scene_;

	TreeMapItem* tree_item_;
	std::vector< VariableItem* > var_items_;

	void UpdateVariableItems(std::vector< CNode* >& selected_nodes, int selected_count, std::vector< QString >& names);

	void UpdateLayout();

	private slots:
	void OnVarSelected(int);
};

#endif