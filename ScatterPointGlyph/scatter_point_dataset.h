#ifndef SCATTER_POINT_DATASET_H_
#define SCATTER_POINT_DATASET_H_

#include <vector>
#include <map>

class ScatterPointDataset
{
public:
	ScatterPointDataset();
	virtual ~ScatterPointDataset();

	void Sample(float left, float right, float bottom, float top);
	void DirectConstruct();
	virtual const char* type() { return "Scatter"; }
	
	int var_num;
	int point_num;
	std::vector< int > sample_index;
	std::map< int, int > node_sample_map;
	std::vector< std::vector< float > > point_pos;
	std::vector< std::vector< float > > point_values;
	std::vector< float > weights;

	std::vector< std::vector< float > > original_point_pos;
	std::vector< std::vector< float > > original_pos_ranges;
	std::vector< std::vector< float > > original_point_values;
	std::vector< std::vector< float > > original_value_ranges;

private:
	void NormalizePosition(std::vector< std::vector< float > >& vec, std::vector< std::vector< float > >& ranges);
	void NormalizeValues(std::vector< std::vector< float > >& vec, std::vector< std::vector< float > >& ranges);
};

#endif