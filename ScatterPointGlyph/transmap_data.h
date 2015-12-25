#ifndef TRANSMAP_DATA_H_
#define TRANSMAP_DATA_H_

#include <vector>
#include "tree_common.h"

class ScatterPointDataset;

class TransMapData
{
public:
	TransMapData();
	~TransMapData();

	// data that must be initialized
	ScatterPointDataset* dataset;
	int cluster_num, var_num;
	std::vector< int > cluster_index;
	std::vector< int > cluster_point_count;
	std::vector< int > cluster_reprentative_color;
	std::vector< int > var_repsentative_color;
	int min_point_num;

	// data that can be completed in ProcessData
	std::vector< CNode* > level_one_nodes;
	std::vector< int > level_one_colors;
	std::vector< CNode* > level_zero_nodes;
	std::vector< int > level_zero_colors;
	std::vector< std::vector< bool > > node_connecting_status;
	std::map< int, CNode* > cluster_node_map;

	// construct the first level nodes and sub-level nodes
	void ClearData();
	void ProcessData();
	void UpdateConnectingStatus();
};

#endif