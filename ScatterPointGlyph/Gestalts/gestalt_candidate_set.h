#ifndef GESTALT_CANDIDATE_SET_H_
#define GESTALT_CANDIDATE_SET_H_

#include <vector>
#include "tree_common.h"

class ScatterPointDataset;

class GestaltCandidateSet
{
public:
	GestaltCandidateSet();
	~GestaltCandidateSet();

	void InitSiteData();
	void ExtractGestaltCandidates(float dis_thresh);
	void ExtractScagnostics();

	// raw point data or sample data
	std::vector< CLeaf* > site_nodes;
	// connecting status of the existing sites
	std::vector< std::vector< bool > > site_connecting_status;
	// existing clusters
	std::vector< CNode* > clusters;
	std::vector< std::vector< bool > > cluster_connecting_status;
	// the cluster index for raw data
	std::vector< int > basic_node_index;

	std::vector< float > var_weights;

	std::vector< std::vector< int > > gestalt_cluster_index;
	std::vector< std::vector< int > >  gestalt_candidates;
	std::vector< std::vector< float > > gestalt_center_pos;
	std::vector< std::vector< float > > gestalt_average_values;

	std::vector< std::vector< double > > gestalt_scagnostics;
};

#endif