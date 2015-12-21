#ifndef TRANSMAP_DATA_H_
#define TRANSMAP_DATA_H_

#include <vector>

class ScatterPointDataset;

class TransMapData
{
public:
	TransMapData();
	~TransMapData();

	ScatterPointDataset* dataset;

	int node_num, var_num;

	std::vector< std::vector< float > > node_center;
	std::vector< std::vector< float > > node_average_value;
	std::vector< std::vector< bool > > node_connecting_status;

	std::vector< int > overall_cluster_index;
	std::vector< int > overall_reprentative_color;

	std::vector< std::vector< int > > var_cluster_index;
	// 3 * var_cluster_index.size()
	std::vector< int > var_repsentative_color;

	void ProcessData();
};

#endif