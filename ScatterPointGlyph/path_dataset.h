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
	std::vector< std::vector< int > > fuzzy_change_values;

	std::vector< bool > is_item_selected;
	std::vector< bool > is_transition_selected;
};

class PathDataset
{
public:
	PathDataset() {}
	~PathDataset() {}

	std::vector< PathRecord* > path_records;
};

#endif