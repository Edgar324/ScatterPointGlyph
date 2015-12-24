#ifndef PATH_DATASET_H_
#define PATH_DATASET_H_

class PathRecord 
{
public:
	PathRecord() {}
	~PathRecord() {}

	std::vector< QColor > item_color;
	std::vector< std::vector< float > > item_values;

	std::vector< std::vector< float > > change_values;

	std::vector< float > var_normalizing_values;
	std::vector< float > var_min_values;
	std::vector< std::string > var_names;
};

class PathDataset
{
public:
	PathDataset() {}
	~PathDataset() {}

	std::vector< PathRecord* > path_records;
};

#endif