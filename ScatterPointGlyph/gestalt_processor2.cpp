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
	: ClusterSolver(), dataset_(NULL), gestalt_candidates_(NULL), linked_tree_(NULL),
	level_num_(2), kmeans_cnum_(600), octree_threshold_(0.3), sampling_mode_(DIRECT), current_level_(0) {
	ProximityExtractor* proximity = new ProximityExtractor;
	SimilarityExtractor* similarity = new SimilarityExtractor;

	property_extractors_.push_back(proximity);
	property_extractors_.push_back(similarity);

	is_property_on_.resize(5, false);
	is_property_on_[0] = true;
	//is_property_on_[1] = true;

	property_thresh_.resize(2);
	property_thresh_[0] = 0.1;
	property_thresh_[1] = 0.3;

	final_label_count_ = 0;
	valid_decreasing_rate_ = 0.3;

	alpha_ = 0.1 / pow(2, level_num_ - 1);
}

GestaltProcessor2::~GestaltProcessor2() {
	if (gestalt_candidates_ != NULL) delete gestalt_candidates_;
}

void GestaltProcessor2::SetData(ScatterPointDataset* data) {
	dataset_ = data;
	if (gestalt_candidates_ != NULL) delete gestalt_candidates_;
	gestalt_candidates_ = new GestaltCandidateSet(dataset_);
	if (linked_tree_ != NULL) delete linked_tree_;
	linked_tree_ = new LinkedTree(dataset_);

	this->SetSamplingMode(sampling_mode_);
}

void GestaltProcessor2::SetSamplingMode(SamplingMode mode) {
	sampling_mode_ = mode;

	switch (sampling_mode_)
	{
	case GestaltProcessor2::DIRECT:
		linked_tree_->ConstructDirectly(level_num_);
		break;
	case GestaltProcessor2::KMEANS:
		linked_tree_->ConstructOnKmeans(level_num_, kmeans_cnum_);
		break;
	case GestaltProcessor2::OCTREE:
		linked_tree_->ConstructOnOctree(level_num_, octree_threshold_);
		break;
	default:
		break;
	}
}

void GestaltProcessor2::SetDisThreshold(float dis_thresh) {
	/*current_level_ = (int)((log(dis_thresh) - log(alpha_)) / log(2.0));
	if (current_level_ >= level_num_) current_level_ = level_num_ - 1;
	if (current_level_ < 0) current_level_ = 0;*/
	current_level_ = 1;
}

void GestaltProcessor2::GetClusterIndex(int& cluster_count, std::vector< int >& cluster_index) {
	cluster_count = linked_tree_->tree_nodes[current_level_].size();
	cluster_index.resize(dataset_->point_pos.size());

	for (int i = 0; i < linked_tree_->tree_nodes[current_level_].size(); ++i) {
		std::queue< Node* > node_queue;
		node_queue.push(linked_tree_->tree_nodes[current_level_][i]);

		while (node_queue.size() != 0) {
			Node* current_node = node_queue.front();
			node_queue.pop();
			if (current_node->type() == Node::LEAF) {
				Leaf* leaf_node = dynamic_cast<Leaf*>(current_node);
				for (int j = 0; j < leaf_node->linked_points.size(); ++j) {
					int point_index = leaf_node->linked_points[j];
					cluster_index[point_index] = i;
				}
			}
			else if (current_node->type() == Node::BRANCH) {
				Branch* branch_node = dynamic_cast<Branch*>(current_node);
				for (int j = 0; j < branch_node->linked_nodes.size(); ++j) node_queue.push(branch_node->linked_nodes[j]);
			}
		}
	}
}

void GestaltProcessor2::UpdateGestaltCandidates(int level) {
	if (linked_tree_->tree_nodes[level].size() != 0 || level <= 0) return;
	if (level > 0 && linked_tree_->tree_nodes[level - 1].size() == 0) return;

	int pre_level = level - 1;

	gestalt_candidates_->InitSiteData(linked_tree_->tree_nodes[pre_level].size());

	for (int i = 0; i < linked_tree_->tree_nodes[pre_level].size(); ++i) {
		std::queue< Node* > node_queue;
		node_queue.push(linked_tree_->tree_nodes[pre_level][i]);

		while (node_queue.size() != 0) {
			Node* current_node = node_queue.front();
			node_queue.pop();
			if (current_node->type() == Node::LEAF) {
				Leaf* leaf_node = dynamic_cast< Leaf* >(current_node);
				gestalt_candidates_->site_point_num[i] += leaf_node->linked_points.size();
				for (int j = 0; j < leaf_node->linked_points.size(); ++j) {
					int point_index = leaf_node->linked_points[j];
					gestalt_candidates_->site_center_pos[i][0] += dataset_->point_pos[point_index][0];
					gestalt_candidates_->site_center_pos[i][1] += dataset_->point_pos[point_index][1];
					gestalt_candidates_->site_average_value[i] += dataset_->point_values[point_index];

					gestalt_candidates_->point_site_id[point_index] = i;
				}
			}
			else if (current_node->type() == Node::BRANCH) {
				Branch* branch_node = dynamic_cast< Branch* >(current_node);
				for (int j = 0; j < branch_node->linked_nodes.size(); ++j) node_queue.push(branch_node->linked_nodes[j]);
			}
		}

		gestalt_candidates_->site_center_pos[i][0] /= gestalt_candidates_->site_point_num[i];
		gestalt_candidates_->site_center_pos[i][1] /= gestalt_candidates_->site_point_num[i];
		gestalt_candidates_->site_average_value[i] /= gestalt_candidates_->site_point_num[i];
	}

	gestalt_candidates_->site_connecting_status = linked_tree_->site_connecting_status;
}

void GestaltProcessor2::run() {
	int pre_level = 0;
	while (pre_level < level_num_ && linked_tree_->tree_nodes[pre_level].size() != 0) pre_level++;
	if (pre_level >= level_num_) return;

	for (int i = pre_level; i <= current_level_; ++i) {
		this->UpdateGestaltCandidates(i);
		this->GenerateCluster(i);
	}
}

void GestaltProcessor2::GenerateCluster(int level) {
	if (gestalt_candidates_ == NULL) return;
	if (linked_tree_->tree_nodes[level].size() != 0) return;

	final_label_.resize(gestalt_candidates_->site_num);
	final_label_.assign(final_label_.size(), -1);

	float dis_thresh = alpha_ * pow(2, level);
	gestalt_candidates_->ExtractGestaltCandidates(dis_thresh);

	labeled_site_count_ = 0;
	final_label_count_ = 0;
	while (labeled_site_count_ < gestalt_candidates_->site_num) {
		ExtractLabels();
		// process valid gestalt
		int property_index, gestalt_index;
		ExtractValidGestalt(property_index, gestalt_index);
		if (property_index == -1  || gestalt_index == -1) break;
		PropertyExtractor* extractor = property_extractors_[property_index];
		for (int i = 0; i < extractor->proposal_gestalt[gestalt_index].size(); ++i) {
			int site_index = extractor->proposal_gestalt[gestalt_index][i];
			final_label_[site_index] = final_label_count_;
			gestalt_candidates_->is_site_labeled[site_index] = true;
		}

		final_label_count_++;
		labeled_site_count_ += extractor->proposal_gestalt[gestalt_index].size();
		std::cout << "Gestalt " << final_label_count_ << "  Point Number: " << extractor->proposal_gestalt[gestalt_index].size() << std::endl;
	}

	for (int i = 0; i < final_label_count_; ++i) {
		Branch* node = new Branch;
		node->seq_index = i;
		linked_tree_->tree_nodes[level].push_back(node);
	}
	for (int i = 0; i < final_label_.size(); ++i) {
		Branch* node = dynamic_cast<Branch*>(linked_tree_->tree_nodes[level][final_label_[i]]);
		node->linked_nodes.push_back(linked_tree_->tree_nodes[level - 1][i]);
	}
	linked_tree_->UpdateConnectingStatus();
}

void GestaltProcessor2::ExtractLabels() {
	try {
		int label_num = gestalt_candidates_->gestalt_candidates.size() + 1 - labeled_site_count_;
		int gestalt_num = gestalt_candidates_->gestalt_candidates.size() - labeled_site_count_;
		int site_num = gestalt_candidates_->site_num; 

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
				extractor->SetData(dataset_, gestalt_candidates_);
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
			}

		delete gc;
	}
	catch (GCException e) {
		e.Report();
	}
}

void GestaltProcessor2::ExtractValidGestalt(int& property_index, int& gestalt_index) {
	std::vector< GestaltInfo > gestalt_infos;

	for (int i = 0; i < is_property_on_.size(); ++i)
		if (is_property_on_[i]) {
			PropertyExtractor* extractor = property_extractors_[i];
			for (int j = 0; j < extractor->proposal_gestalt.size(); ++j) {
				GestaltInfo info_item;
				info_item.property_index = i;
				info_item.gestalt_index = j;

				info_item.values[0] = 0;
				for (int k = 0; k < extractor->proposal_gestalt[j].size(); ++k) {
					int site_index = extractor->proposal_gestalt[j][k];
					info_item.values[0] += extractor->data_cost[site_index][j];
				}

				info_item.values[1] = 0;
				for (int k = 0; k < extractor->proposal_gestalt[j].size(); ++k)
					if (extractor->result_label[k] == j) info_item.values[1] += 1;
				info_item.values[1] = (extractor->proposal_gestalt[j].size() - info_item.values[1]) / extractor->proposal_gestalt[j].size();

				gestalt_infos.push_back(info_item);
			}
		}
	
	SortGestaltInfos(gestalt_infos, 0, gestalt_infos.size() - 1, 1);
	int end_index = 1;
	while (end_index < gestalt_infos.size() && gestalt_infos[end_index].values[1] < valid_decreasing_rate_) end_index++;
	SortGestaltInfos(gestalt_infos, 0, end_index - 1, 0);

	property_index = gestalt_infos[0].property_index;
	gestalt_index = gestalt_infos[0].gestalt_index;
}

void GestaltProcessor2::SortGestaltInfos(std::vector< GestaltInfo >& infos, int begin, int end, int sort_index) {
	if (begin >= end) return;

	int first = begin, last = end;
	GestaltInfo key = infos[first];
	while (first < last) {
		while (first < last && infos[last].values[sort_index] >= key.values[sort_index]) --last;
		infos[first] = infos[last];
		while (first < last && infos[first].values[sort_index] <= key.values[sort_index]) ++first;
		infos[last] = infos[first];
	}
	infos[first] = key;
	SortGestaltInfos(infos, begin, first - 1, sort_index);
	SortGestaltInfos(infos, first + 1, end, sort_index);
}