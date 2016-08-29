#ifndef SCATTER_POINT_DATASET_H_
#define SCATTER_POINT_DATASET_H_

#include <vector>
#include <map>
using namespace std;
#include <QtCore/QString>
#include <QtGui/QColor>

class ScatterPointDataset
{
public:
	ScatterPointDataset();
	virtual ~ScatterPointDataset();

    enum DataType {
        POINT_DATA = 0x0,
        GRID_DATA,
        VOLUME_DATA
    };

    virtual DataType type() { return POINT_DATA; }

	// original data, var_names, var_weights, original_point_values should be set manually
	int var_num = {0};
	vector<QString> var_names;
	vector<float> var_weights;
    vector<QColor> var_colors;

    int point_num = {0};

	vector<vector<float>> original_point_values;
    // the value ranges of the points before normalization
	vector<vector<float>> original_value_ranges;

	// mds result
	vector<vector<float>> original_point_pos;
	// the position ranges of the points before normalization
	vector<vector<float>> original_pos_ranges;
    float max_pos_range;

    // normalized point positions
	vector<vector<float>> normalized_point_pos;
	// normalized point variable values
	vector<vector<float>> normalized_point_values;

    vector<bool> is_valid;

	
	void ApplyMds();
    void ApplyTsne();
    void ApplyNormal(int axis_one, int axis_two);

	// Direct use the normalized original point positions and variable values
	// TODO: Provide sampling methods for large scale data
	void DirectConstruct();
    void ConstructWithValueRanges(vector<vector<float>>& ranges);
    void NormalizePos();
	void ClearData();

protected:
    vector<int> selected_vars;
	// The x and y coordinate share the same scale. All values are normalized to [0, 1]
	void NormalizePosition(vector<vector<float>>& vec, vector<vector<float>>& ranges);
	void NormalizeValues(vector<vector<float>>& vec, vector<vector<float>>& ranges);
};

#endif