#ifndef TOUR_PATH_GENERATOR_H_
#define TOUR_PATH_GENERATOR_H_

#include <vector>

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

	PathRecord* GetPath();

	std::vector< int > edge_list;
	std::vector< int > tour_node_list;

private:
	TransMapData* trans_data_;
};

#endif