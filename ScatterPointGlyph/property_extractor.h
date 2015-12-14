#ifndef PROPERTY_EXTRACTOR_H_
#define PROPERTY_EXTRACTOR_H_

#include <vector>

#include "gestalt_candidate_set.h"

class PropertyExtractor
{
public:
	PropertyExtractor();
	~PropertyExtractor();

	void SetData(GestaltCandidateSet* candidates);
	virtual void ExtractCosts(float thres);
	void ExtractFitness();
	void SetParameters(float data_cost_scale, float data_cost_bias, float smooth_cost_scale, float smooth_cost_bias, float label_cost_scale, float label_cost_bias);

	std::vector< float > label_cost;
	std::vector< std::vector< float > > smooth_cost;
	std::vector< std::vector< float > > data_cost;

	std::vector< std::vector< int > > proposal_gestalt;
	std::vector< std::vector< int > > proposal_clusters;
	std::vector< float > fitness;
	std::vector< int > result_label;

protected:
	GestaltCandidateSet* gestalt_candidates;
	const float maximum_cost = { 100.0 };
	float scales[3];
	float bias[3];

	virtual void ExtractProposalGestalt(float thres);
	void NormalizeVec(std::vector< float >& vec);
	void NormalizeVec(std::vector< std::vector< float > >& vec);
};

#endif