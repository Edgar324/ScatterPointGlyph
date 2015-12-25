#include "tour_path_generator.h"
#include "path_dataset.h"
#include "transmap_data.h"

TourPathGenerator::TourPathGenerator() {

}

TourPathGenerator::~TourPathGenerator() {

}

void TourPathGenerator::SetData(TransMapData* data) {
	trans_data_ = data;
}

bool TourPathGenerator::GenerateRoundPath() {
	return true;
}

bool TourPathGenerator::GeneratePath(int begin_index /* = -1 */, int end_index /* = -1 */) {
	if (begin_index == -1 || end_index == -1) return false;
	return true;
}

PathRecord* TourPathGenerator::GetPath() {
	if (trans_data_ == 0) return 0;

	PathRecord* record = new PathRecord;
	record->item_values.resize(trans_data_->level_one_nodes.size());
	for (int i = 0; i < trans_data_->level_one_nodes.size(); ++i) {
		record->item_values[i] = trans_data_->level_one_nodes[i]->average_values;
	}
	record->change_values.resize(trans_data_->level_one_nodes.size() - 1);
	for (int i = 0; i < trans_data_->level_one_nodes.size() - 1; ++i) {
		CNode* current_node = trans_data_->level_one_nodes[i];
		CNode* next_node = trans_data_->level_one_nodes[i + 1];
		record->change_values[i].resize(current_node->average_values.size());
		for (int j = 0; j < current_node->average_values.size(); ++j)
			record->change_values[i][j] = next_node->average_values[j] - current_node->average_values[j];
	}
	for (int k = 0; k < trans_data_->level_one_nodes[0]->average_values.size(); ++k) 
		record->var_names.push_back("Test");

	return record;
}