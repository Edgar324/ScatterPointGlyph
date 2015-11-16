#ifndef PROXIMITY_EXTRACTOR_H_
#define PROXIMITY_EXTRACTOR_H_

#include "property_extractor.h"

class ProximityExtractor : public PropertyExtractor
{
public:
	ProximityExtractor();
	~ProximityExtractor();

	void ExtractCosts(float thres);

private:
	void ExtractProposalGestalt(float thres);
};

#endif