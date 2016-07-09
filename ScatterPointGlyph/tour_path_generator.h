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

	static bool GenerateRoundPath(std::vector<CNode*>& nodes, std::vector<int>& tour_list);
	static bool GenerateRoundPath(std::vector<std::vector<float>>& node_dis, std::vector<int>& tour_list);
};

#endif