#ifndef TOUR_PATH_GENERATOR_H_
#define TOUR_PATH_GENERATOR_H_

class TransMapData;
class PathRecord;

class TourPathGenerator
{
public:
	TourPathGenerator();
	~TourPathGenerator();

	void SetData(TransMapData* data);
	bool GenerateRoundPath();
	bool GeneratePath(int begin_index = -1, int end_index = -1);

	PathRecord* GetPath();

private:
	TransMapData* trans_data_;
};

#endif