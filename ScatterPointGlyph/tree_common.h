#ifndef TREE_COMMON_H_
#define TREE_COMMON_H_

#include <QThread>
#include "cnode.h"

class ScatterPointDataset;

class TreeCommon : public QThread
{
public:
	TreeCommon(ScatterPointDataset* data);
	virtual ~TreeCommon();

	enum TreeMode 
	{
		EXPLORATION_MODE = 0x0,
		VIEWING_MODE
	};

	CBranch* root() { return root_; }
	ScatterPointDataset* data() { return dataset_; }
	void SetTreeMode(TreeMode mode);

	// Functions for generating the leaf nodes
	void ConstructDirectly();

	void GetClusterResult(int level, int& cluster_num, std::vector<int>& cluster_index);
	void GetClusterResult(int level, std::vector<CNode*>& level_nodes);
	int GetMaxLevel() { return this->max_level_; }

	// Sort the tree nodes according to the node_ids
	void SortTree(std::vector<int>& node_ids);

	// Tree node operations
	void SplitCluster(int cluster_index);
	void MergeClusters(std::vector<int>& cluster_index);

	void GetNodeValues(CNode* node, int var_index, std::vector<float>& values);

	void Traverse(CNode* node, std::vector<int>& linked_points);
	void Traverse(int level, std::vector<CNode*>& nodes);
	void TraversAllNodes(CNode* root_node, std::vector<CNode*>& nodes);
	void TraverseLevelNodes(int level, CNode* root_node, std::vector<CNode*>& nodes);

protected:
	ScatterPointDataset* dataset_;
	CBranch* root_;
	int max_level_;
	float min_edge_length_;
	float data_dis_scale_;
	TreeMode tree_mode_;

	std::map< int, CNode*> id_node_map_;
	CNode* common_parent_node_;

    std::map< int, int > grid_id_seq_map_;

    std::vector<std::vector<bool>> node_connecting_status_;

	void run();
	virtual void GenerateClusters() = 0;
	virtual void BeginClustering() = 0;
	virtual void SplitNode(CBranch* node) = 0;

    void GetConnectionStatus(std::vector<CNode*>& nodes, std::vector<std::vector<bool>>& connecting_status, float& min_edge_length);

	// Update the level of the node and its children
	void ResetLevel(CNode* node, int level);

	void AssignLeafLevel(CNode* node, int level);

	// Assign Tree Color
	// This algorithm is based on the paper: 
	// Tennekes, Martijn, and Edwin de Jonge. "Tree colors: color schemes for tree-structured data." 
	// Visualization and Computer Graphics, IEEE Transactions on 20, no. 12 (2014): 2072-2081.
	void AssignColor(CNode* node, float hstart, float hend, float factor = 0.75, bool perm = true, bool rev = true);

	// Assign sequence index for each branch node
	void ResetSortingIndex(CNode* node);
	int SortNode(CNode* node, std::vector<int>& node_ids, int& node_count);

	// Tree node operations
	// Update node parameters
	void ProgressNodeData(CNode* node);
	void ProgressNodeAndParentData(CNode* node);

	void RemoveChildNode(CNode* node, bool is_empty_deleted = false);
	void UpdateChildLevel(CBranch* node);
	int FindCommonParent(CNode* node, std::vector<int>& node_ids);

    // TODO: Add other construction methods
	void ConstructOnSlic(float thre);
};

#endif