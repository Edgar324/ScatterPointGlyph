#ifndef GESTALT_CANDIDATE_SET_H_
#define GESTALT_CANDIDATE_SET_H_

#include <vector>

class ScatterPointDataset;

class GestaltCandidateSet
{
public:
	GestaltCandidateSet(ScatterPointDataset* data);
	~GestaltCandidateSet();

	void ExtractGestaltCandidates();

	std::vector< int > point_site_id;
	std::vector< int > site_point_num;
	std::vector< std::vector< float > > site_center_pos;
	std::vector< float > site_average_value;
	std::vector< bool > is_site_labeled;

	std::vector< std::vector< int > >  gestalt_candidates;
	std::vector< std::vector< bool > > connecting_status;
	std::vector< std::vector< int > > edge_links;

private:
	ScatterPointDataset* dataset_;
};

#endif