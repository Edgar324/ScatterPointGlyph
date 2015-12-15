#ifndef IMMEDIATE_TREE_H_
#define IMMEDIATE_TREE_H_

#include "tree_common.h"

class ImmediateTree : public TreeCommon
{
public:
	ImmediateTree(ScatterPointDataset* data);
	~ImmediateTree();

	void SetRadiusThreshold(float max_radius);
	void SetSimilarityThreshold(float thre);
	void SetProximityThreshold(float thre);
	void SetSampleSize(int num);

	virtual void GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index);
	virtual void GetClusterResult(float dis_per_piexl, int& cluster_num, std::vector< int >& cluster_index);

protected:
	virtual void run();
	virtual void GenerateCluster();

	float similarity_threshold_, proximity_threshold_;
	float max_radius_threshold_;
	bool is_threshold_updated_;
	int sample_num_;

	std::vector< int > node_cluster_index_;
};

#endif