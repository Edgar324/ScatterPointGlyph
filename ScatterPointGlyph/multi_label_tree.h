#ifndef MULTI_LABEL_TREE_H_
#define MULTI_LABEL_TREE_H_

#include "tree_common.h"

class MultiLabelProcessor;

class MultiLabelTree : public TreeCommon
{
public:
	MultiLabelTree(ScatterPointDataset* data);
	~MultiLabelTree();

	// Set radius of level 1
	void SetRadiusThreshold(float max_radius);
	void SetUncertaintyThreshold(float un_threshold);
	int GetRadiusLevel(float radius);

protected:
	virtual void GenerateClusters();
	virtual void SplitNode(CBranch* node);

	// Generate new graph based on the cluster results based on each variable
	void GenerateSegmentUncertainty(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status, std::vector< std::vector< float > >& edge_weight);

	MultiLabelProcessor* processor_;

	float max_radius_threshold_;
	float un_threshold_;
	float factor_ = 1.5;
};

#endif