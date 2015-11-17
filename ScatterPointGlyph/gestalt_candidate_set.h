#ifndef GESTALT_CANDIDATE_SET_H_
#define GESTALT_CANDIDATE_SET_H_

#include <vector>

class ScatterPointDataset;

class GestaltCandidateSet
{
public:
	GestaltCandidateSet(ScatterPointDataset* data);
	~GestaltCandidateSet();

	void ExtractGestaltCandidates(float dis_thresh);

	std::vector< int > point_site_id;
	int site_num;
	std::vector< int > site_point_num;
	std::vector< std::vector< float > > site_center_pos;
	std::vector< float > site_average_value;
	std::vector< std::vector< bool > > site_connecting_status;
	std::vector< bool > is_site_labeled;

	std::vector< std::vector< int > >  gestalt_candidates;

private:
	ScatterPointDataset* dataset_;

	void InitSiteData();
	void Triangulation();
};

#endif