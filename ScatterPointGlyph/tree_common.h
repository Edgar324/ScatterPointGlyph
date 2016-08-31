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

    enum TreeType {
        HIERARCHICAL_TREE = 0x0,
        NCUTS_TREE,
        MULTI_LABEL_TREE,
        VIEW_DEPENDENT_TREE,
        CLUSTER_PROJECTION_TREE
    };

	CBranch* root() { return root_; }
	ScatterPointDataset* data() { return point_dataset_; }

    void GetNodes(float left, float right, float bottom, float top, float glyph_radius,
        vector<CNode*>& pre_level_nodes, vector<CNode*>& current_level_nodes);
    void GetNodes(int level, vector<CNode*>& level_nodes);
    CNode* GetNode(int node_id);
    void GetNodePoints(int node_id, vector<int>& point_ids);
    void GetNodePoints(CNode* node, std::vector<int>& point_ids);
    int GetMaxLevel() { return this->max_level_; }

    void MergeNodes(vector<int>& node_ids);
    void SplitNodeOnce(int node_id);
    void SplitNodeRecursively(int node_id, float std_dev_threshold);

    void Clear();

    virtual void AutoConstructTree(float std_dev_threshold) = 0;
    virtual void SplitNode(CBranch* node) = 0;
    virtual TreeType type() = 0;
    virtual void ConstructTree(float left, float right, float bottom, float top, float glyph_radius) {}

protected:
    ScatterPointDataset* point_dataset_ = NULL;
	CBranch* root_ = NULL;
	int max_level_ = -1;
    float data_dis_scale_ = 0.5;

    const int MAX_ALLOWED_LEVEL = 7;

    map<int, CNode*> id_node_map_;

    void run();

    // Tree node operations
    void ConstructRootNode();
    void ConstructLargeScaleRootNode();
	// Update node parameters
	void ProgressNodeData(CNode* node);

   	// Update the level of the node and its children
	void ResetLevel(CNode* node, int level);


//public:
//    void GetNodePoints(int level, int& cluster_num, vector<int>& cluster_index);
//	// Functions for generating the leaf nodes
//
//	// Sort the tree nodes according to the node_ids
//	void SortTree(std::vector<int>& node_ids);
//
//	// Tree node operations
//	void SplitCluster(int cluster_index);
//	void MergeClusters(std::vector<int>& cluster_index);
//
//	void GetNodeValues(CNode* node, int var_index, std::vector<float>& values);
//
//	
//	void Traverse(int level, std::vector<CNode*>& nodes);
//	void TraversAllNodes(CNode* root_node, std::vector<CNode*>& nodes);
//	void TraverseLevelNodes(int level, CNode* root_node, std::vector<CNode*>& nodes);
//
//protected:
//
//	float min_edge_length_;
//
//	
//	CNode* common_parent_node_;
//
//    std::map< int, int > grid_id_seq_map_;
//
//    std::vector<std::vector<bool>> node_connecting_status_;
//
//	virtual void GenerateClusters() = 0;
//	virtual void BeginClustering() = 0;
//
//    void GetConnectionStatus(std::vector<CNode*>& nodes, std::vector<std::vector<bool>>& connecting_status, float& min_edge_length);
//

//
//	void AssignLeafLevel(CNode* node, int level);
//
//	// Assign Tree Color
//	// This algorithm is based on the paper: 
//	// Tennekes, Martijn, and Edwin de Jonge. "Tree colors: color schemes for tree-structured data." 
//	// Visualization and Computer Graphics, IEEE Transactions on 20, no. 12 (2014): 2072-2081.
//	void AssignColor(CNode* node, float hstart, float hend, float factor = 0.75, bool perm = true, bool rev = true);
//
//	// Assign sequence index for each branch node
//	void ResetSortingIndex(CNode* node);
//	int SortNode(CNode* node, std::vector<int>& node_ids, int& node_count);
//
//	void RemoveChildNode(CNode* node, bool is_empty_deleted = false);
//	void UpdateChildLevel(CBranch* node);
//	int FindCommonParent(CNode* node, std::vector<int>& node_ids);
//
//    // TODO: Add other construction methods
//	void ConstructOnSlic(float thre);
};

#endif