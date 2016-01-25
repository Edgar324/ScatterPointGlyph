#include "continuity_extractor.h"

#include <fstream>
#include <iostream>

#define PRE_LOAD

ContinuityExtractor::ContinuityExtractor() {

}

ContinuityExtractor::~ContinuityExtractor() {

}

void ContinuityExtractor::SetData(std::vector< float >& x, std::vector< float >& y) {
	xPos = x;
	yPos = y;

	ExtractContinuity();
}

void ContinuityExtractor::ExtractContinuity() {
	float max_dis = -1e10;
	for (int i = 0; i < xPos.size() - 1; ++i)
		for (int j = i + 1; j < xPos.size(); ++j) {
			float temp_dis = sqrt(pow(xPos[i] - xPos[j], 2) + pow(yPos[i] - yPos[j], 2));
			if (temp_dis > max_dis) max_dis = temp_dis;
		}

#ifndef PRE_LOAD
	continuity.resize(xPos.size());
	for (int i = 0; i < xPos.size(); ++i) {
		continuity[i].resize(xPos.size());
		memset(continuity[i].data(), 0, sizeof(float) * xPos.size());
		continuity[i][i] = 1e05;
	}

	int step_count = 5;
	float current_scale = 1;
	std::vector< std::vector< bool > > is_connected;
	is_connected.resize(xPos.size());
	for (int t = 0; t < step_count; ++t) {
		current_scale += pow(2, t);
		float current_dis = max_dis * pow(0.5, t) * 0.1;
		for (int i = 0; i < is_connected.size(); ++i){
			is_connected[i].resize(xPos.size());
			is_connected[i].assign(xPos.size(), false);
			is_connected[i][i] = true;
		}
		for (int i = 0; i < xPos.size() - 1; ++i)
			for (int j = i + 1; j < xPos.size(); ++j) {
				float temp_dis = sqrt(pow(xPos[i] - xPos[j], 2) + pow(yPos[i] - yPos[j], 2));
				if (temp_dis < current_dis) is_connected[i][j] = true;
				is_connected[j][i] = is_connected[i][j];
			}
		for (int i = 0; i < xPos.size(); ++i){
			std::cout << i << std::endl;
			for (int j = 0; j < xPos.size(); ++j)
				for (int k = 0; k < xPos.size(); ++k)
					if (!is_connected[j][k] && is_connected[j][i] && is_connected[i][k]) {
						is_connected[j][k] = true;
						is_connected[k][j] = true;
						break;
					}
		}
		for (int i = 0; i < xPos.size() - 1; ++i)
			for (int j = i + 1; j < xPos.size(); ++j) {
				if (is_connected[i][j]) continuity[i][j] += current_scale;
				continuity[j][i] = continuity[i][j]; 
			}
	}

	float max_scale = -1e10;
	for (int i = 0; i < xPos.size() - 1; ++i)
		for (int j = i + 1; j < xPos.size(); ++j) {
			if (i != j && continuity[i][j] > max_scale) max_scale = continuity[i][j];
		}
	for (int i = 0; i < xPos.size() - 1; ++i)
		for (int j = i + 1; j < xPos.size(); ++j){
			continuity[i][j] /= max_scale;
			continuity[j][i] = continuity[i][j];
		}
	for (int i = 0; i < xPos.size(); ++i) continuity[i][i] = 1;

	std::ofstream output_stream("./TestData/test.txt", std::ios::binary);
	for (int i = 0; i < xPos.size(); ++i)
		output_stream.write((char*)continuity[i].data(), continuity[i].size() * sizeof(float));
	output_stream.close();
#else
	continuity.resize(xPos.size());
	for (int i = 0; i < xPos.size(); ++i) {
		continuity[i].resize(xPos.size());
	}
	std::ifstream input_stream("./TestData/test.txt", std::ios::binary);
	for (int i = 0; i < xPos.size(); ++i)
		input_stream.read((char*)continuity[i].data(), xPos.size() * sizeof(float));
	input_stream.close();
#endif // !PRE_LOAD

}