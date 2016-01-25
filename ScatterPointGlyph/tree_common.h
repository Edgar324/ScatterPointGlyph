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
	void ConstructOnOctree(float thre);
	void ConstructOnKmeans(int basic_cnum);
	void ConstructOnRandomSample(int sample_num);
	void ConstructDirectly();

	void GetClusterResult(int level, int& cluster_num, std::vector< int >& cluster_index);
	void GetClusterResult(int level, std::vector< CNode* >& level_nodes);

	// Sort the tree nodes according to the node_ids
	void SortTree(std::vector< int >& node_ids);

protected:
	ScatterPointDataset* dataset_;

	CBranch* root_;
	float min_edge_length_;
	int max_level_;

	virtual void run();
	virtual void GenerateCluster(CBranch* node = NULL);

	void InitializeSortingIndex();
	int SortNode(CNode* node, std::vector< int >& node_ids, int& node_count);
	
	void Traverse(CNode* node, std::vector< int >& linked_points);
	void Traverse(int level, std::vector< CNode* >& nodes);
	void Traverse(int level, CNode* root, std::vector< CNode* >& nodes);

	void ProgressNode(CNode* node);
	void ResetLevel(CNode* node, int level);

	void AssignColor(CNode* node, float hstart, float hend, float factor = 0.75, bool perm = true, bool rev = true);
};

#endif