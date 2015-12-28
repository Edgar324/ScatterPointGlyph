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

private:
	TransMapData* trans_data_;

	std::vector< int > edge_list_;
	std::vector< int > tour_node_list_;
};

#endif