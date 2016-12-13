#ifndef TREE_COMMON_H_
#define TREE_COMMON_H_

#include <QThread>
#include "cnode.h"

class MultivariateDataset;

class TreeCommon : public QThread
{
public:
	TreeCommon(MultivariateDataset* dataset);
	virtual ~TreeCommon();

    enum TreeType {
        GEO_HIERARCHICAL_TREE = 0x0,
        GEO_NCUTS_TREE,
        GEO_MULTI_LABEL_TREE,
        GEO_VIEW_DEPENDENT_MULTI_LABEL_TREE,
        KNN_MULTI_LABEL_TREE,
        KNN_VIEW_DEPENDENT_MULTI_LABEL_TREE
    };

	CBranch* root() { return root_; }
	MultivariateDataset* mv_dataset() { return mv_dataset_; }
    int max_level() { return this->max_level_; }

    // Return nodes that are in a specified level. Hierarchical control.
    void GetNodes(int level, vector<CNode*>& level_nodes);

    // Return nodes that are best fit for the view. View dependent control.
    void GetNodes(float left, float right, float bottom, float top, vector<CNode*>& nodes);

    // Return a node according to its id
    CNode* GetNode(int node_id);

    virtual TreeType type() = 0;

    // Split a branch node using specified algorithm.
    // If the node has branch children, all children will be destroyed.
    // Only nodes with all leaf children are available for splitting
    virtual void SplitNode(CBranch* node) = 0;

    // Construct a hierarchical tree based on the std_dev_threshold
    virtual void AutoConstructTree(float std_dev_threshold = 0.0) = 0;

    // Construct a two-level tree based on the view dependent control
    virtual void ConstructTree(float left, float right, float bottom, float top, float glyph_radius) = 0;

    // Clear the tree and generate the basic two-level tree with root and leaf nodes
    virtual void Clear() = 0;

protected:
    MultivariateDataset* mv_dataset_ = NULL;

	CBranch* root_ = NULL;
    vector<CLeaf*> leaf_nodes_;
    map<int, CNode*> id_node_map_;

	int max_level_ = 0;
    const int MAX_ALLOWED_LEVEL = 7;
    
    virtual void BuildLeafs() = 0;

    void Initialize();
    void run();


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