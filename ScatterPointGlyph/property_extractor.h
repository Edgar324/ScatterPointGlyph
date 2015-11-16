#ifndef PROPERTY_EXTRACTOR_H_
#define PROPERTY_EXTRACTOR_H_

#include <vector>

#include "gestalt_candidate_set.h"
#include "scatter_point_dataset.h"

class PropertyExtractor
{
public:
	PropertyExtractor();
	~PropertyExtractor();

	void SetData(ScatterPointDataset* data, GestaltCandidateSet* candidates);
	virtual void ExtractCosts(float thres);

	std::vector< float > label_cost;
	std::vector< std::vector< float > > smooth_cost;
	std::vector< std::vector< float > > data_cost;

	std::vector< std::vector< int > > proposal_gestalt;

	std::vector< int > result_label;

protected:
	ScatterPointDataset* dataset;
	GestaltCandidateSet* gestalt_candidates;

	virtual void ExtractProposalGestalt(float thres);
	void NormalizeVec(std::vector< float >& vec);
	void NormalizeVec(std::vector< std::vector< float > >& vec);
};

#endif