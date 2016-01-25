#ifndef SIMILARITY_EXTRACTOR_H_
#define SIMILARITY_EXTRACTOR_H_

#include "property_extractor.h"

class SimilarityExtractor : public PropertyExtractor
{
public:
	SimilarityExtractor();
	~SimilarityExtractor();

	virtual void ExtractCosts(float thres);

private:
	void ExtractProposalGestalt(float thres);
};

#endif