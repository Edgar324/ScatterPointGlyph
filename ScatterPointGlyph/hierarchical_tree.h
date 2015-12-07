#ifndef HIERARCHICAL_TREE_H_
#define HIERARCHICAL_TREE_H_

#include "tree_common.h"

class ScatterPointDataset;

class HierarchicalTree : public TreeCommon
{
public:
	HierarchicalTree(ScatterPointDataset* data);
	~HierarchicalTree();

	void GenerateCluster(float dis_per_pixel = 0.01, int min_pixel_radius = 1);
	void SetVariableWeights(std::vector< float >& weights);

private:
	int max_level_;
	float similarity_threshold_, proximity_threshold_;
	float fitness_threshold_;
	std::vector< int > node_cluster_index_;

	std::vector< float > var_weights_;
};

#endif