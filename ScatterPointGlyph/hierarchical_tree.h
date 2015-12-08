#ifndef HIERARCHICAL_TREE_H_
#define HIERARCHICAL_TREE_H_

#include "tree_common.h"

class ScatterPointDataset;

class HierarchicalTree : public TreeCommon
{
public:
	HierarchicalTree(ScatterPointDataset* data);
	~HierarchicalTree();

	void SetMaxLevel(int level);
	void SetSimilarityThreshold(float thre);
	void SetProximityThreshold(float thre);
	void SetFitnessThreshold(float thre);
	void SetRadiusParameters(float min_pixel_radius);

	virtual void GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index);
	virtual void GetClusterResult(float dis_per_piexl, int& cluster_num, std::vector< int >& cluster_index);

protected:
	virtual void run();

private:
	int max_level_;
	float similarity_threshold_, proximity_threshold_;
	float fitness_threshold_;
	float min_pixel_radius_;

	std::vector< int > node_cluster_index_;

	void GenerateCluster(int min_pixel_radius);
};

#endif