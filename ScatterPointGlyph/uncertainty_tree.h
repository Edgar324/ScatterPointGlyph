#ifndef UNCERTAINTY_TREE_H_
#define UNCERTAINTY_TREE_H_

#include "tree_common.h"

class MultiLabelProcessor;

class UncertaintyTree : public TreeCommon
{
public:
	UncertaintyTree(ScatterPointDataset* data);
	~UncertaintyTree();

	void SetRadiusThreshold(float max_radius);
	void SetSampleSize(int num);

	virtual void GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index);
	virtual void GetClusterResult(float dis_per_piexl, int& cluster_num, std::vector< int >& cluster_index);
	std::vector< float >& GetUncertainty();

protected:
	virtual void run();
	virtual void GenerateCluster();
	void GenerateUncertainty();

	MultiLabelProcessor* processor_;

	float max_radius_threshold_;
	bool is_threshold_updated_;
	int sample_num_;

	std::vector< std::vector< float > > edge_weights_;
	std::vector< float > leaf_node_un_;
	std::vector< float > un_value_;
};

#endif