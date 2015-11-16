#ifndef SCATTER_POINT_DATASET_H_
#define SCATTER_POINT_DATASET_H_

#include <vector>

class ScatterPointDataset
{
public:
	ScatterPointDataset() {}
	~ScatterPointDataset() {}

	std::vector< std::vector< float > > point_pos;
	std::vector< float > point_values;

	std::vector< std::vector< float > > origin_point_values;
	std::vector< float > weights;
};

#endif