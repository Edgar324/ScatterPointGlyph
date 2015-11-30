#ifndef GESTALT_CANDIDATE_SET_H_
#define GESTALT_CANDIDATE_SET_H_

#include <vector>

class ScatterPointDataset;

class GestaltCandidateSet
{
public:
	GestaltCandidateSet(ScatterPointDataset* data);
	~GestaltCandidateSet();

	void InitSiteData(int num);
	void ExtractGestaltCandidates(float dis_thresh);
	void ExtractScagnostics();

	std::vector< int > point_site_id;

	int site_num;
	std::vector< int > site_point_num;
	std::vector< std::vector< float > > site_center_pos;
	std::vector< float > site_average_value;
	std::vector< bool > is_site_labeled;
	std::vector< std::vector< bool > > site_connecting_status;

	std::vector< std::vector< int > >  gestalt_candidates;
	std::vector< std::vector< double > > gestalt_scagnostics;

private:
	ScatterPointDataset* dataset_;
};

#endif