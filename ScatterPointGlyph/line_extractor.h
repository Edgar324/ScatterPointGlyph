#ifndef LINE_EXTRACTOR_H_
#define LINE_EXTRACTOR_H_

#include <vector>

class LineExtractor
{
public:
	LineExtractor();
	~LineExtractor();

	void ExtractLines(std::vector< double >& x, std::vector< double >& y, std::vector< float >& paras);
};

#endif