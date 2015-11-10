#include "line_extractor.h"

#include <iostream>
#include <mrpt/math/ransac_applications.h>
using namespace mrpt;
using namespace mrpt::math;

LineExtractor::LineExtractor() {

}

LineExtractor::~LineExtractor() {

}

void LineExtractor::ExtractLines(std::vector< double >& x, std::vector< double >& y, std::vector< float >& paras) {
	std::vector< std::pair< size_t, TLine2D > > detectedLines;
	const double DIST_THRESHOLD = 200;

	CVectorDouble xs, ys;
	xs.resize(x.size());
	ys.resize(y.size());
	for (int i = 0; i < x.size(); ++i) {
		xs[i] = x[i];
		ys[i] = y[i];
	}
	ransac_detect_2D_lines(xs, ys, detectedLines, DIST_THRESHOLD, 20);
	std::cout << "RANSAC method: ransac_detect_2D_lines" << std::endl;
	std::cout << " " << detectedLines.size() << " lines detected." << std::endl;
	for (int i = 0; i < detectedLines.size(); ++i) {
		paras.push_back(detectedLines[i].second.coefs[0]);
		paras.push_back(detectedLines[i].second.coefs[1]);
		paras.push_back(detectedLines[i].second.coefs[2]);
	}
}