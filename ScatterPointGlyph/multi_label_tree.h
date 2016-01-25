#ifndef MULTI_LABEL_TREE_H_
#define MULTI_LABEL_TREE_H_

#include "tree_common.h"

class MultiLabelProcessor;

class MultiLabelTree : public TreeCommon
{
public:
	MultiLabelTree(ScatterPointDataset* data);
	~MultiLabelTree();

	void SetRadiusThreshold(float max_radius);
	void SetUncertaintyThreshold(float un_threshold);
	void SetSampleSize(int num);
	
	void SetSegmentUncertaintyOn();
	void SetSegmentUncertaintyOff();

	void SplitCluster(int cluster_index);
	void MergeClusters(std::vector< int >& cluster_index);
	void AddUserDefinedCluster(int origin_cluster, std::vector< int >& point_index);

protected:
	virtual void GenerateClusters();
	virtual void GenerateCluster(CBranch* node);
	void GenerateSegmentUncertainty(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status, std::vector< std::vector< float > >& edge_weight);

	void ProgressNodeAndParent(CNode* node);
	void RemoveChildNode(CNode* node, bool is_empty_deleted = false);
	void UpdateChildLevel(CBranch* node);
	int FindCommonParent(CNode* node, std::vector< int >& node_ids);

	MultiLabelProcessor* processor_;
	std::map< int, CNode* > id_node_map_;

	CNode* common_parent_node_;

	float max_radius_threshold_;
	float un_threshold_;
	int sample_num_;
	float factor_ = 1.5;

	bool is_segment_uncertainty_applied_;
};

#endif