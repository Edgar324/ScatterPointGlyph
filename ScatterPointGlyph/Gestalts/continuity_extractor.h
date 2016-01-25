#ifndef CONTINUITY_EXTRACTOR_H_
#define CONTINUITY_EXTRACTOR_H_

#include <vector>

class ContinuityExtractor
{
public:
	ContinuityExtractor();
	~ContinuityExtractor();

	void SetData(std::vector< float >& x, std::vector< float >& y);

	std::vector< std::vector< float > > continuity;

private:
	std::vector< float > xPos;
	std::vector< float > yPos;

	void ExtractContinuity();
};

#endif