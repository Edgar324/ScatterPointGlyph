#ifndef GRID_SCATTER_POINT_DATASET_H_
#define GRID_SCATTER_POINT_DATASET_H_

#include "scatter_point_dataset.h"

class GridScatterPointDataset : public ScatterPointDataset
{
public:
	GridScatterPointDataset();
	virtual ~GridScatterPointDataset();

	virtual const char* type() { return "Grid"; }

	int w, h;
};

#endif