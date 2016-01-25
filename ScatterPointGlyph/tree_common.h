#ifndef TREE_COMMON_H_
#define TREE_COMMON_H_

#include <QThread>
#include "cnode.h"

class ScatterPointDataset;

//#define RADAR_GLYPH

class TreeCommon : public QThread
{
public:
	TreeCommon(ScatterPointDataset* data);
	virtual ~TreeCommon();

	CBranch* root() { return root_; }

	// Functions for generating the leaf nodes
	void ConstructDirectly();
	void ConstructOnKmeans(int basic_cnum);
	// TODO: Add other construction methods
	void ConstructOnOctree(float thre);
	void ConstructOnRandomSample(int sample_num);

	void GetClusterResult(int level, int& cluster_num, std::vector< int >& cluster_index);
	void GetClusterResult(int level, std::vector< CNode* >& level_nodes);

	// Sort the tree nodes according to the node_ids
	void SortTree(std::vector< int >& node_ids);

protected:
	ScatterPointDataset* dataset_;
	CBranch* root_;
	float min_edge_length_;
	float data_dis_scale_;

	void run();
	virtual void GenerateClusters() = 0;

	// Update node parameters
	void ProgressNode(CNode* node);

	// Update the level of the node and its children
	void ResetLevel(CNode* node, int level);

	// Assign Tree Color
	// This algorithm is based on the paper: 
	// Tennekes, Martijn, and Edwin de Jonge. "Tree colors: color schemes for tree-structured data." 
	// Visualization and Computer Graphics, IEEE Transactions on 20, no. 12 (2014): 2072-2081.
	void AssignColor(CNode* node, float hstart, float hend, float factor = 0.75, bool perm = true, bool rev = true);

	void InitializeSortingIndex();
	int SortNode(CNode* node, std::vector< int >& node_ids, int& node_count);

	void Traverse(CNode* node, std::vector< int >& linked_points);
	void Traverse(int level, std::vector< CNode* >& nodes);
	void Traverse(int level, CNode* root, std::vector< CNode* >& nodes);
};

#endif