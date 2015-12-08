#include "property_extractor.h"

 
PropertyExtractor::PropertyExtractor() {

}

PropertyExtractor::~PropertyExtractor() {

}

void PropertyExtractor::SetData(GestaltCandidateSet* candidates) {
	gestalt_candidates = candidates;
}

void PropertyExtractor::ExtractCosts(float thres) {

}

void PropertyExtractor::ExtractProposalGestalt(float thres) {

}

void PropertyExtractor::ExtractFitness() {
	fitness.resize(this->proposal_gestalt.size());
	fitness.assign(this->proposal_gestalt.size(), 0);

	std::vector< float > label_fitness;
	label_fitness.resize(this->proposal_clusters.size());

	for (int i = 0; i < this->proposal_gestalt.size(); ++i) {
		label_fitness.assign(label_fitness.size(), 0);
		for (int j = 0; j < this->proposal_gestalt[i].size(); ++j) {
			int site_index = this->proposal_gestalt[i][j];
			label_fitness[result_label[gestalt_candidates->basic_node_index[site_index]]] += 1;
		}
		for (int j = 0; j < label_fitness.size(); ++j) {
			label_fitness[j] /= this->proposal_gestalt[i].size();
			if (label_fitness[j] > fitness[i]) fitness[i] = label_fitness[j];
		}
	}
}

void PropertyExtractor::NormalizeVec(std::vector< float >& vec) {
	float minValue = 1e10;
	float maxValue = -1e10;

	for (int i = 0; i < vec.size(); ++i){
		if (minValue > vec[i] && vec[i] >= 0) minValue = vec[i];
		if (maxValue < vec[i] && vec[i] >= 0) maxValue = vec[i];
	}

	if (maxValue - minValue != 0 && maxValue - minValue < 1e10) {
		for (int i = 0; i < vec.size(); ++i)
			if (vec[i] >= 0)
				vec[i] = vec[i] / maxValue;
			else
				vec[i] = 10.0;
	}
	else {
		for (int i = 0; i < vec.size(); ++i) vec[i] = 1.0;
	}
}

void PropertyExtractor::NormalizeVec(std::vector< std::vector< float > >& vec) {
	float minValue = 1e10;
	float maxValue = -1e10;

	for (int i = 0; i < vec.size(); ++i)
		for (int j = 0; j < vec[i].size(); ++j) {
			if (minValue > vec[i][j] && vec[i][j] >= 0) minValue = vec[i][j];
			if (maxValue < vec[i][j] && vec[i][j] >= 0) maxValue = vec[i][j];
		}

	if (maxValue - minValue != 0 && maxValue - minValue < 1e10) {
		for (int i = 0; i < vec.size(); ++i)
			for (int j = 0; j < vec[i].size(); ++j)
				if (vec[i][j] >= 0)
					vec[i][j] = vec[i][j] / maxValue;
				else
					vec[i][j] = 10.0;
	}
	else {
		for (int i = 0; i < vec.size(); ++i)
			for (int j = 0; j < vec[i].size(); ++j)
				vec[i][j] = 1.0;
	}
}