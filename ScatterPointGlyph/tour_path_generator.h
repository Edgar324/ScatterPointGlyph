#ifndef TOUR_PATH_GENERATOR_H_
#define TOUR_PATH_GENERATOR_H_

#include <vector>
#include "tree_common.h"

class TransMapData;
class PathRecord;

class TourPathGenerator
{
public:
	TourPathGenerator();
	~TourPathGenerator();

	void SetData(TransMapData* data);
	bool GenerateRoundPath();
	bool GenerateSpanningTree();
	bool GeneratePath(int begin_index = -1, int end_index = -1);
	
	static bool GenerateRoundPath(std::vector< CNode* >& nodes, std::vector< int >& tour_list);
	static bool GenerateRoundPath(std::vector< std::vector< float > >& node_dis, std::vector< int >& tour_list);

	PathRecord* GetPath();

	std::vector< int > edge_list;
	std::vector< int > tour_node_list;

private:
	TransMapData* trans_data_;
};

#endif