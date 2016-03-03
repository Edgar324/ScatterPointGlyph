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

    enum DataType {
        POINT_DATA = 0x0,
        GRID_DATA
    };

    virtual DataType type() { return POINT_DATA; }

	// original data, var_names, var_weights, original_point_values should be set manually
	int var_num = {0};
	std::vector< QString > var_names;
	std::vector< float > var_weights;
	std::vector< std::vector< float > > original_point_values;

	int point_num = {0};
	// normalized point positions
	std::vector< std::vector< float > > point_pos;
	// normalized point variable values
	std::vector< std::vector< float > > point_values;

	// mds result
	std::vector< std::vector< float > > original_point_pos;
	// the position ranges of the points before normalization
	std::vector< std::vector< float > > original_pos_ranges;
    float max_pos_range;

	// the value ranges of the points before normalization
	std::vector< std::vector< float > > original_value_ranges;

	void ManualSelectDim(std::vector< bool >& is_dim_selected);
	// TODO: Add automatic dimension reduction method
	void AutoDimReduction(int dim_num);

	// Execute MDS on the normalized multi-variate data
	// The implementation of MDS is from Quan Wang.
	// More information about this function can be retrieved from "SimpleMatrix.h" or
	// https://sites.google.com/site/simpmatrix/
	void ExecMds();

	// Direct use the normalized original point positions and variable values
	// TODO: Provide sampling methods for large scale data
	void DirectConstruct();
	void ClearData();

protected:
	// The x and y coordinate share the same scale. All values are normalized to [0, 1]
	void NormalizePosition(std::vector< std::vector< float > >& vec, std::vector< std::vector< float > >& ranges);
	void NormalizeValues(std::vector< std::vector< float > >& vec, std::vector< std::vector< float > >& ranges);
};

#endif