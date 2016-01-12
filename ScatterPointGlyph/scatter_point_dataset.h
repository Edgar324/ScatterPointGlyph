#ifndef SCATTER_POINT_DATASET_H_
#define SCATTER_POINT_DATASET_H_

#include <vector>
#include <map>
#include <QtCore/QString>

class ScatterPointDataset
{
public:
	ScatterPointDataset();
	virtual ~ScatterPointDataset();

	void Sample(float left, float right, float bottom, float top);
	void DirectConstruct();
	virtual const char* type() { return "Scatter"; }
	
	// sampled data
	int var_num;
	int point_num;
	std::vector< int > sample_index;
	std::vector< std::vector< float > > point_pos;
	std::vector< std::vector< float > > point_values;

	// mds result
	std::vector< std::vector< float > > original_point_pos;
	std::vector< std::vector< float > > original_pos_ranges;

	// original data
	std::vector< QString > var_names;
	std::vector< float > var_weights;
	std::vector< std::vector< float > > original_point_values;
	std::vector< std::vector< float > > original_value_ranges;

	void AutoDimReduction(int dim_num);
	void SelectDim(std::vector< bool >& is_dim_selected);

	void ExecMds();

private:
	void NormalizePosition(std::vector< std::vector< float > >& vec, std::vector< std::vector< float > >& ranges);
	void NormalizeValues(std::vector< std::vector< float > >& vec, std::vector< std::vector< float > >& ranges);
};

#endif