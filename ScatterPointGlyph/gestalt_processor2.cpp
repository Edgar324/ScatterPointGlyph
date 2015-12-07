#include "gestalt_processor2.h"
#include <iostream>
#include <queue>
#include <time.h>

#include "scatter_point_dataset.h"
#include "gestalt_candidate_set.h"
#include "similarity_extractor.h"
#include "proximity_extractor.h"
#include "linked_tree.h"
#include "./gco-v3.0/GCoptimization.h"
#include "./gco-v3.0/energy.h"

GestaltProcessor2::GestaltProcessor2() 
	: ClusterSolver() {
	ProximityExtractor* proximity = new ProximityExtractor;
	SimilarityExtractor* similarity = new SimilarityExtractor;

	property_extractors_.push_back(proximity);
	property_extractors_.push_back(similarity);

	is_property_on_.resize(5, false);
	is_property_on_[0] = true;
	is_property_on_[1] = true;

	property_thresh_.resize(2);
	property_thresh_[0] = 0.1;
	property_thresh_[1] = 0.3;

	gestalt_candidates_ = new GestaltCandidateSet;

	final_label_count_ = 0;
}

GestaltProcessor2::~GestaltProcessor2() {
	if (gestalt_candidates_ != NULL) delete gestalt_candidates_;
}

void GestaltProcessor2::SetData(std::vector< CLeaf* >& nodes, std::vector< CNode* >& clusters, 
	std::vector< int >& node_index, std::vector< std::vector< bool > >& node_connecting_status,
	std::vector< float >& var_weights) {
	gestalt_candidates_->site_nodes = nodes;
	gestalt_candidates_->clusters = clusters;
	gestalt_candidates_->basic_node_index = node_index;
	gestalt_candidates_->site_connecting_status = node_connecting_status;
	gestalt_candidates_->var_weights = var_weights;
	gestalt_candidates_->InitSiteData();
}

void GestaltProcessor2::SetDisThreshold(float dis_thresh) {
	dis_threshold_ = dis_thresh;
}

void GestaltProcessor2::SetPropertyOn(GestaltProperty property) {
	is_property_on_[property] = true;
}

void GestaltProcessor2::SetPropertyOff(GestaltProperty property) {
	is_property_on_[property] = false;
}

void GestaltProcessor2::SetThreshold(GestaltProperty property, float thresh) {
	property_thresh_[property] = thresh;
}

void GestaltProcessor2::GetCluster(int& cluster_count, std::vector< int >& cluster_index) {
	
}

void GestaltProcessor2::GetCluster(float fitness, std::vector< int >& cluster_index) {
	float best_fitness = -1e10;
	int property_index = -1;
	int gestalt_index = -1;
	for (int i = 0; i < property_extractors_.size(); ++i)
		if (is_property_on_[i]) {
			PropertyExtractor* extractor = property_extractors_[i];
			for (int j = 0; j < extractor->fitness.size(); ++j)
				if (extractor->fitness[j] > fitness && extractor->fitness[j] > best_fitness && extractor->proposal_clusters[j].size() >= 2) {
					best_fitness = extractor->fitness[j];
					property_index = i;
					gestalt_index = j;
				}
		}

	cluster_index.clear();
	if (property_index != -1 && gestalt_index != -1) {
		PropertyExtractor* extractor = property_extractors_[property_index];
		for (int i = 0; i < extractor->proposal_clusters[gestalt_index].size(); ++i)
			cluster_index.push_back(extractor->proposal_clusters[gestalt_index][i]);
	}
}

void GestaltProcessor2::run() {
	this->GenerateCluster();
}

void GestaltProcessor2::GenerateCluster() {
	if (gestalt_candidates_ == NULL) return;
	
	gestalt_candidates_->ExtractGestaltCandidates(dis_threshold_);
	this->ExtractLabels();
}

void GestaltProcessor2::ExtractLabels() {
	try {
		int label_num = gestalt_candidates_->gestalt_candidates.size();
		int site_num = gestalt_candidates_->clusters.size(); 

		// Step 1: construct class
		GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(site_num, label_num);

		// Step 2: multi-label graph cut
		for (int i = 0; i < site_num - 1; ++i)
			for (int j = i + 1; j < site_num; ++j) {
				if (gestalt_candidates_->site_connecting_status[i][j]) gc->setNeighbors(i, j);
			}

		for (int p = 0; p < is_property_on_.size(); ++p)
			if (is_property_on_[p]) {
				PropertyExtractor* extractor = property_extractors_[p];
				extractor->SetData(gestalt_candidates_);
				extractor->ExtractCosts(property_thresh_[p]);

				for (int i = 0; i < site_num; ++i) {
					for (int j = 0; j < label_num; ++j) {
						gc->setDataCost(i, j, (int)(extractor->data_cost[i][j] * 10000));
					}
				}
				for (int i = 0; i < label_num; ++i)
					for (int j = 0; j < label_num; ++j)
						gc->setSmoothCost(i, j, (int)(extractor->smooth_cost[i][j] * 10000));

				std::vector< int > temp_label_cost;
				temp_label_cost.resize(label_num);
				for (int i = 0; i < label_num; ++i) temp_label_cost[i] = (int)(extractor->label_cost[i] * 10000);
				gc->setLabelCost(temp_label_cost.data());

				//std::cout << "Property: " << p << "    Before optimization energy is " << gc->compute_energy() << std::endl;
				gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
				//std::cout << "Property: " << p << "    After optimization energy is " << gc->compute_energy() << std::endl;

				extractor->result_label.resize(site_num);
				for (int i = 0; i < site_num; ++i) extractor->result_label[i] = gc->whatLabel(i);

				extractor->ExtractFitness();
			}

		delete gc;
	}
	catch (GCException e) {
		e.Report();
	}
}