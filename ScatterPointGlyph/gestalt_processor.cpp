#include "gestalt_processor.h"
#include <time.h>
#include <queue>
#include <iostream>
#include <vtkDelaunay2D.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkLine.h>
#include <vtkCellLocator.h>
#include <vtkSmartPointer.h>
#include <vtkDelaunay2D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>
#include "./gco-v3.0/GCoptimization.h"

GestaltProcessor::GestaltProcessor() 
	: cluster_num_(100), dis_scale_(0.5), distance_threshold_(0.1){
	is_sample_needed_ = false;
}

GestaltProcessor::~GestaltProcessor() {

}

void GestaltProcessor::SetData(std::vector< std::vector< float > >& pos, std::vector< std::vector< float > >& values, std::vector< float >& weight){
	point_pos_ = pos;
	point_values_.resize(pos.size());
	for (int i = 0; i < pos.size(); ++i) {
		point_values_[i] = 0;
		for (int j = 0; j < weight.size(); ++j)
			point_values_[i] += values[i][j] * weight[j];
	}
}

void GestaltProcessor::SetGestaltThreshold(std::vector< float >& threshold){
	gestalt_threshold_ = threshold;
}

void GestaltProcessor::SetDistanceThreshold(float thresh) {
	distance_threshold_ = thresh;
}

void GestaltProcessor::GetFinalGlyphPoint(std::vector< int >& point_index) {
	point_index.clear();
	for (int i = 0; i < point_pos_.size(); ++i)
		if (final_label_[cluster_index_[i]] == final_gestalt_count_) point_index.push_back(i);

	std::cout << "Glyph point number: " << point_index.size() << std::endl;
}

void GestaltProcessor::GetKmeansClusterResult(std::vector< std::vector< int > >& index){
	index.clear();
	index.resize(cluster_num_);
	for (int i = 0; i < point_pos_.size(); ++i)
		index[cluster_index_[i]].push_back(i);
}

void GestaltProcessor::GenerateLayout() {
	// step one: kmeans to generate a series of candidate sites
	KmeansCluster();
	emit KMeansClusterFinished();

	// step two: generate connecting graph using delaunay triangulation
	Triangulation();

	is_node_labeled_.resize(cluster_num_);
	is_node_labeled_.assign(cluster_num_, false);

	final_gestalt_count_ = -1;
	final_label_.resize(cluster_num_);
	final_label_.assign(cluster_num_, -1);

	int rule_number = 2;
	gestalt_info_.resize(rule_number * cluster_num_);
	for (int i = 0; i < gestalt_info_.size(); ++i) gestalt_info_[i].resize(4);
	remaining_gestalt_node_map_.resize(rule_number);
	remaining_node_gestalt_map_.resize(rule_number);

	proposal_gestalt_.resize(rule_number);
	optimized_labels_.resize(rule_number);

	bool is_all_labeld = false;
	while (!is_all_labeld) {
		gestalt_info_.clear();

		// step three: generate gestalts based on different rules
		// 1) Similarity   2) Proximity  3) Continuity  4) closure
		// step four: generate optimization result based on the nodes
		ExtractNodeDistance();
		ProcessGestaltRule(0);
		ProcessGestaltRule(1);

		bool is_one_run_passed = false;
		// step five: remove the strongest gestalt, goto step four until all nodes are removed
		float valid_decrease_rate = 0.3;
		bool is_gestalt_valid;
		int begin_index = 0;
		do {
			Sort(gestalt_info_, begin_index, gestalt_info_.size() - 1, 1);
			int end_index = begin_index;
			while (end_index < gestalt_info_.size() && gestalt_info_[end_index][1] < valid_decrease_rate) end_index++;
			is_gestalt_valid = true;
			if (end_index > begin_index) Sort(gestalt_info_, begin_index, end_index, 0);
			else is_gestalt_valid = false;

			if (is_one_run_passed && !is_gestalt_valid) break;

			// pick the best gestalt
			final_gestalt_count_++;
			int selected_gestalt_rule = (int)gestalt_info_[begin_index][2];
			int gestalt_index = (int)gestalt_info_[begin_index][3];
			for (int i = 0; i < proposal_gestalt_[selected_gestalt_rule][gestalt_index].size(); ++i) {
				int node_index = proposal_gestalt_[selected_gestalt_rule][gestalt_index][i];
				if (is_node_labeled_[node_index]) continue;
				is_node_labeled_[node_index] = true;
				final_label_[proposal_gestalt_[selected_gestalt_rule][gestalt_index][i]] = final_gestalt_count_;
			}

			emit FinalGlyphUpdated();

			// update the decreasing rate and check whether valid gestalt exist
			for (int rule = 0; rule < rule_number; ++rule) {
				std::vector< float > temp_decrease_rate;
				temp_decrease_rate.resize(proposal_gestalt_[rule].size(), 0);
				for (int i = 0; i < proposal_gestalt_[rule].size(); ++i){
					temp_decrease_rate[i] = 0;
					for (int j = 0; j < proposal_gestalt_[rule][i].size(); ++j){
						int index = proposal_gestalt_[rule][i][j];
						if (!is_node_labeled_[index] && optimized_labels_[rule][remaining_node_gestalt_map_[rule][index]] == i) temp_decrease_rate[i] += 1;
					}
				}

				for (int i = 0; i < proposal_gestalt_[rule].size(); ++i){
					temp_decrease_rate[i] = (proposal_gestalt_[rule][i].size() - temp_decrease_rate[i]) / proposal_gestalt_[rule][i].size();
				}
				for (int i = begin_index + 1; i < gestalt_info_.size(); ++i) {
					int temp_rule_index = (int)gestalt_info_[i][2];
					int temp_gestalt_index = (int)gestalt_info_[i][3];
					if (rule == temp_rule_index) gestalt_info_[i][1] = temp_decrease_rate[temp_gestalt_index];
				}
			}

			begin_index++;
			is_one_run_passed = true;
		} while (is_gestalt_valid && begin_index < gestalt_info_.size());

		int labeled_count = 0;
		for (int i = 0; i < is_node_labeled_.size(); ++i)
			if (is_node_labeled_[i]) labeled_count++;
		if (cluster_num_ - labeled_count <= 2) {
			for (int i = 0; i < is_node_labeled_.size(); ++i)
				if (!is_node_labeled_[i]) {
					final_gestalt_count_++;
					final_label_[i] = final_gestalt_count_;
				}
			is_all_labeld = true;
		}
	}
}

void GestaltProcessor::KmeansCluster() {
	if (!is_sample_needed_) {
		cluster_num_ = point_pos_.size();
		cluster_center_.resize(cluster_num_);
		cluster_node_count_.resize(cluster_num_, 1);
		cluster_index_.resize(point_pos_.size());
		for (int i = 0; i < cluster_num_; ++i) {
			cluster_center_[i].resize(3);
			cluster_center_[i][0] = point_pos_[i][0];
			cluster_center_[i][1] = point_pos_[i][1];
			cluster_center_[i][2] = point_values_[i];
		}
		for (int i = 0; i < point_pos_.size(); ++i) cluster_index_[i] = i;
	}
	else {
		cluster_num_ = point_pos_.size() / 3;
		/// Generate random center
		cluster_center_.resize(cluster_num_);
		cluster_node_count_.resize(cluster_num_);

		std::vector< bool > is_selected;
		is_selected.resize(point_values_.size(), false);

		/// remove duplicated data
		for (int i = 0; i < point_pos_.size() - 1; ++i)
			if (!is_selected[i]) {
				for (int j = i + 1; j < point_pos_.size(); ++j) {
					float temp_dis = abs(point_pos_[i][0] - point_pos_[j][0]) + abs(point_pos_[i][1] - point_pos_[j][1]);
					if (temp_dis <= 1e-20) is_selected[j] = true;
				}
			}
		srand((unsigned int)time(0));
		for (int i = 0; i < cluster_num_; ++i) {
			int temp_index;
			do {
				temp_index = (int)((float)rand() / RAND_MAX * (point_values_.size() - 1) + 0.499);
			} while (is_selected[temp_index]);

			cluster_center_[i].resize(3);
			cluster_center_[i][0] = point_pos_[temp_index][0];
			cluster_center_[i][1] = point_pos_[temp_index][1];
			cluster_center_[i][2] = point_values_[temp_index];
			is_selected[temp_index] = true;
		}

		cluster_index_.resize(point_pos_.size());
		for (int i = 0; i < cluster_index_.size(); ++i) cluster_index_[i] = -1;

		/// Update cluster
		bool is_center_updated;
		int iteration_count = 0;
		do {
			is_center_updated = false;
			for (int i = 0; i < point_pos_.size(); ++i) {
				float min_dis_index = -1;
				float min_dis = 1e20;
				for (int j = 0; j < cluster_num_; ++j) {
					float temp_dis = 0;
					temp_dis += sqrt(pow(point_pos_[i][0] - cluster_center_[j][0], 2) + pow(point_pos_[i][1] - cluster_center_[j][1], 2)) * dis_scale_;
					temp_dis += abs(point_values_[j] - point_values_[i]) * (1.0 - dis_scale_);
					if (temp_dis < min_dis) {
						min_dis = temp_dis;
						min_dis_index = j;
					}
				}
				if (cluster_index_[i] != min_dis_index) {
					is_center_updated = true;
					cluster_index_[i] = min_dis_index;
				}
			}
			if (!is_center_updated) break;

			/// update cluster center
			for (int i = 0; i < cluster_num_; ++i)
				memset(cluster_center_[i].data(), 0, cluster_center_[i].size() * sizeof(float));
			memset(cluster_node_count_.data(), 0, cluster_node_count_.size() * sizeof(int));
			for (int i = 0; i < cluster_index_.size(); ++i) {
				cluster_node_count_[cluster_index_[i]]++;
				cluster_center_[cluster_index_[i]][0] += point_pos_[i][0];
				cluster_center_[cluster_index_[i]][1] += point_pos_[i][1];
				cluster_center_[cluster_index_[i]][2] += point_values_[i];
			}
			for (int i = 0; i < cluster_num_; ++i)
				if (cluster_node_count_[i] != 0) {
					for (int j = 0; j < cluster_center_[i].size(); ++j)
						cluster_center_[i][j] /= cluster_node_count_[i];
				}

			iteration_count++;
		} while (is_center_updated && iteration_count < 30);

		for (int i = 0; i < cluster_num_; ++i) {
			if (cluster_node_count_[i] == 0) {
				while (cluster_node_count_[i] == 0){
					for (int j = i; j < cluster_num_ - 1; ++j){
						cluster_center_[j] = cluster_center_[j + 1];
						cluster_node_count_[j] = cluster_node_count_[j + 1];
					}
					for (int j = 0; j < cluster_index_.size(); ++j)
						if (cluster_index_[j] > i) cluster_index_[j] -= 1;
					cluster_num_--;
				}
			}
		}
		cluster_center_.resize(cluster_num_);
		cluster_node_count_.resize(cluster_num_);
	}
	
	
}

void GestaltProcessor::Triangulation() {
	vtkSmartPointer< vtkPoints > points = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i < cluster_num_; ++i) 
		points->InsertNextPoint(cluster_center_[i][0], cluster_center_[i][1], 0);
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(polydata);
	delaunay->Update();
	vtkPolyData* triangle_out = delaunay->GetOutput();

	connect_status_.resize(cluster_num_);
	for (int i = 0; i < cluster_num_; ++i) connect_status_[i].resize(cluster_num_, false);
	vtkIdTypeArray* idarray = triangle_out->GetPolys()->GetData();
	for (int i = 0; i < triangle_out->GetNumberOfPolys(); ++i){
		double* temp = idarray->GetTuple(i * 3);
		int id1 = (int)(temp[0]);
		temp = idarray->GetTuple(i * 3 + 1);
		int id2 = (int)(temp[0]);
		temp = idarray->GetTuple(i * 3 + 2);
		int id3 = (int)(temp[0]);
		connect_status_[id1][id2] = true;
		connect_status_[id2][id1] = true;
		connect_status_[id2][id3] = true;
		connect_status_[id3][id2] = true;
		connect_status_[id1][id3] = true;
		connect_status_[id3][id1] = true;
	}
}

void GestaltProcessor::ProcessGestaltRule(int rule) {
	// step one: extract proposal gestalt
	// in this method, we use the cluster center as a representative gestalt center,
	// all nodes can be added to the gestalt according to the gestalt rule
	ExtractProposalGestalt(rule);

	// step two: label optimization
	// step three: generate gestalt energy
	ExecLabelOptimization(rule);
}

void GestaltProcessor::ExtractNodeDistance() {
	node_distance_.resize(cluster_num_);
	for (int i = 0; i < cluster_num_; ++i) node_distance_[i].resize(cluster_num_);

	std::vector< bool > is_reached;
	is_reached.resize(cluster_num_);

	for (int i = 0; i < cluster_num_; ++i) {
		if (is_node_labeled_[i]) continue;
		node_distance_[i].assign(cluster_num_, 1e10);
		is_reached.assign(cluster_num_, false);
		node_distance_[i][i] = 0;

		for (int j = 0; j < is_node_labeled_.size(); ++j)
			if (is_node_labeled_[j]) {
				is_reached[j] = true;
			}

		for (int j = 0; j < cluster_num_ - 1; ++j) {
			float min_dist = 1e20;
			int min_index = -1;
			for (int k = 0; k < cluster_num_; ++k)
				if (!is_reached[k] && node_distance_[i][k] < min_dist) {
					min_dist = node_distance_[i][k];
					min_index = k;
				}
			if (min_index == -1) break;
			is_reached[min_index] = true;

			for (int k = 0; k < cluster_num_; ++k)
				if (!is_reached[k] && connect_status_[min_index][k]) {
					float temp_dis = sqrt(pow(cluster_center_[j][0] - cluster_center_[k][0], 2) + pow(cluster_center_[j][1] - cluster_center_[k][1], 2));
					if (node_distance_[i][k] > node_distance_[i][min_index] + temp_dis) node_distance_[i][k] = node_distance_[i][min_index] + temp_dis;
				}
		}
	}
}

void GestaltProcessor::UpdateDecreasingRate(int rule) {

}

void GestaltProcessor::ExtractProposalGestalt(int rule){
	proposal_gestalt_[rule].clear();

	std::vector< bool > is_reached;
	is_reached.resize(cluster_num_);
	remaining_node_gestalt_map_[rule].clear();
	remaining_gestalt_node_map_[rule].clear();

	int remaining_node_count = 0;
	for (int i = 0; i < cluster_num_; ++i) {
		if (is_node_labeled_[i]) continue;
		remaining_node_gestalt_map_[rule].insert(std::map< int, int >::value_type(i, remaining_node_count));
		remaining_gestalt_node_map_[rule].insert(std::map< int, int >::value_type(remaining_node_count, i));
		remaining_node_count++;

		std::vector< int > gestalt;
		for (int j = 0; j < node_distance_[i].size(); ++j)
			if (node_distance_[i][j] < distance_threshold_) {
				float rule_value = 0;
				switch (rule)
				{
				case 0:
					rule_value = abs(cluster_center_[j][2] - cluster_center_[i][2]);
					break;
				case 1:
					rule_value = sqrt(pow(cluster_center_[j][0] - cluster_center_[i][0], 2) + pow(cluster_center_[j][1] - cluster_center_[i][1], 2));
					break;
				default:
					break;
				}
				if (rule_value < gestalt_threshold_[rule]) gestalt.push_back(j);
			}

		/// TODO: filter out identical gestalt
		proposal_gestalt_[rule].push_back(gestalt);
	}
}

void GestaltProcessor::ExecLabelOptimization(int rule) {
	try {
		int gestalt_num = proposal_gestalt_[rule].size();
		// Step 1: construct class
		GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(gestalt_num, gestalt_num);

		// Step 2: construct data cost, smooth cost and label cost for each label
		std::vector< std::vector< float > > data_cost, smooth_cost;
		std::vector< float > label_cost;
		data_cost.resize(gestalt_num);
		smooth_cost.resize(gestalt_num);
		label_cost.resize(gestalt_num, -1);
		for (int i = 0; i < gestalt_num; ++i) {
			data_cost[i].resize(gestalt_num, -1);
			smooth_cost[i].resize(gestalt_num, -1);
		}

		// label cost
		std::vector< std::vector< float > > gestalt_average;
		std::vector< int > gestalt_node_count;
		std::vector< float > gestalt_variance;
		gestalt_average.resize(gestalt_num);
		gestalt_node_count.resize(gestalt_num, 0);
		gestalt_variance.resize(gestalt_num, 0);
		gestalt_average.resize(gestalt_num);
		int average_length = 1;
		if (rule == 1) average_length = 2;
		for (int i = 0; i < gestalt_num; ++i) gestalt_average[i].resize(average_length, 0);

		for (int i = 0; i < proposal_gestalt_[rule].size(); ++i) {
			for (int j = 0; j < proposal_gestalt_[rule][i].size(); ++j){
				int index = proposal_gestalt_[rule][i][j];
				gestalt_node_count[i] += cluster_node_count_[index];
				switch (rule){
				case 0:
					gestalt_average[i][0] += cluster_center_[index][2] * cluster_node_count_[index];
					break;
				case 1:
					gestalt_average[i][0] += cluster_center_[index][0] * cluster_node_count_[index];
					gestalt_average[i][1] += cluster_center_[index][1] * cluster_node_count_[index];
					break;
				default:
					break;
				}
			}
		}
		for (int i = 0; i < gestalt_num; ++i)
			for (int j = 0; j < average_length; ++j)
				gestalt_average[i][j] /= gestalt_node_count[i];

		for (int i = 0; i < gestalt_num; ++i) {
			for (int j = 0; j < proposal_gestalt_[rule][i].size(); ++j){
				for (int k = 0; k < point_pos_.size(); ++k)
					if (cluster_index_[k] == proposal_gestalt_[rule][i][j]) {
						switch (rule){
						case 0:
							gestalt_variance[i] += pow(point_values_[k] - gestalt_average[i][0], 2);
							break;
						case 1:
							gestalt_variance[i] += pow(point_pos_[k][0] - gestalt_average[i][0], 2) + pow(point_pos_[k][1] - gestalt_average[i][1], 2);
							break;
						default:
							break;
						}
					}
			}
		}
		for (int i = 0; i < gestalt_num; ++i)
			gestalt_variance[i] = sqrt(gestalt_variance[i] / gestalt_node_count[i]);
		for (int i = 0; i < label_cost.size(); ++i) label_cost[i] = gestalt_variance[i];

		// data cost
		for (int i = 0; i < proposal_gestalt_[rule].size(); ++i) {
			for (int j = 0; j < proposal_gestalt_[rule][i].size(); ++j) {
				int temp_index = proposal_gestalt_[rule][i][j];
				int gestalt_index = remaining_node_gestalt_map_[rule][temp_index];
				switch (rule){
				case 0:
					data_cost[gestalt_index][i] = abs(cluster_center_[temp_index][2] - gestalt_average[i][0]);
					break;
				case 1:
					data_cost[gestalt_index][i] = sqrt(pow(cluster_center_[temp_index][0] - gestalt_average[i][0], 2) + pow(cluster_center_[temp_index][1] - gestalt_average[i][1], 2));
					break;
				default:
					break;
				}
			}
		}
		for (int i = 0; i < gestalt_num; ++i){
			data_cost[i][i] = 0;
		}

		// smooth cost
		for (int i = 0; i < gestalt_num; ++i) {
			for (int j = i + 1; j < gestalt_num; ++j) {
				int i_index = remaining_gestalt_node_map_[rule][i];
				int j_index = remaining_gestalt_node_map_[rule][j];
				float temp_dis = sqrt(pow(cluster_center_[i_index][0] - cluster_center_[j_index][0], 2) + pow(cluster_center_[i_index][1] - cluster_center_[j_index][1], 2));
				if (temp_dis < distance_threshold_) {
					smooth_cost[i][j] = 0.8;
				}
				else
					smooth_cost[i][j] = 0.4;
				//smooth_cost[i][j] = 1.0 / (temp_dis + 1);
				smooth_cost[j][i] = smooth_cost[i][j];
				
			}
			smooth_cost[i][i] = 0;
		}

		//if (!CheckMetric(smooth_cost)) exit(-1);

		NormalizeVec(data_cost);
		NormalizeVec(label_cost);
		//NormalizeVec(smooth_cost);

		// Step 3: multi-label graph cut
		for (int i = 0; i < gestalt_num - 1; ++i)
			for (int j = i + 1; j < gestalt_num; ++j) {
				int i_index = remaining_gestalt_node_map_[rule][i];
				int j_index = remaining_gestalt_node_map_[rule][j];
				if (connect_status_[i_index][j_index]) gc->setNeighbors(i, j);
			}
		for (int i = 0; i < gestalt_num; ++i) {
			for (int j = 0; j < gestalt_num; ++j) {
				gc->setDataCost(i, j, data_cost[i][j]);
				gc->setSmoothCost(i, j, smooth_cost[i][j]);
			}
			gc->setLabelCost(label_cost[i]);
		}

		std::cout << "Before optimization energy is " << gc->compute_energy() << std::endl;
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		std::cout << "After optimization energy is " << gc->compute_energy() << std::endl;

		optimized_labels_[rule].resize(gestalt_num);
		for (int i = 0; i < gestalt_num; i++)
			optimized_labels_[rule][i] = gc->whatLabel(i);

		int pre_count = gestalt_info_.size();
		gestalt_info_.resize(pre_count + gestalt_num);
		for (int i = pre_count; i < pre_count + gestalt_num; ++i) gestalt_info_[i].resize(4);
		// update energy
		for (int i = 0; i < proposal_gestalt_[rule].size(); ++i){
			gestalt_info_[pre_count + i][0] = 0;
			for (int j = 0; j < proposal_gestalt_[rule][i].size(); ++j){
				int index = remaining_node_gestalt_map_[rule][proposal_gestalt_[rule][i][j]];
				gestalt_info_[pre_count + i][0] += data_cost[index][i];
			}
		}

		// update gestalt decrease rate
		for (int i = 0; i < proposal_gestalt_[rule].size(); ++i){
			gestalt_info_[pre_count + i][1] = 0;
			for (int j = 0; j < proposal_gestalt_[rule][i].size(); ++j){
				int index = remaining_node_gestalt_map_[rule][proposal_gestalt_[rule][i][j]];
				if (optimized_labels_[rule][index] == i) gestalt_info_[pre_count + i][1] += 1;
			}
		}

		for (int i = 0; i < proposal_gestalt_[rule].size(); ++i){
			gestalt_info_[pre_count + i][1] = (proposal_gestalt_[rule][i].size() - gestalt_info_[pre_count + i][1]) / proposal_gestalt_[rule][i].size();
		}

		for (int i = 0; i < proposal_gestalt_[rule].size(); ++i) {
			gestalt_info_[pre_count + i][2] = rule;
			gestalt_info_[pre_count + i][3] = i;
		}

		delete gc;
	}
	catch (GCException e) {
		e.Report();
	}
}

void GestaltProcessor::Sort(std::vector< std::vector< float > >& value, int begin, int end, int sort_index) {
	std::vector< float > temp_vec;
	for (int i = begin; i < end - 1; ++i)
		for (int j = i + 1; j < end; ++j)
			if (value[i][sort_index] > value[j][sort_index]) {
				temp_vec = value[i];
				value[i] = value[j];
				value[j] = temp_vec;
			}
}

void GestaltProcessor::NormalizeVec(std::vector< float >& vec) {
	float minValue = 1e10;
	float maxValue = -1e10;

	for (int i = 0; i < vec.size(); ++i){
		if (minValue > vec[i] && vec[i] >= 0) minValue = vec[i];
		if (maxValue < vec[i] && vec[i] >= 0) maxValue = vec[i];
	}

	if (maxValue - minValue != 0 && maxValue - minValue < 1e10) {
		for (int i = 0; i < vec.size(); ++i)
			if (vec[i] >= 0)
				vec[i] = (vec[i] - minValue) / (maxValue - minValue);
			else
				vec[i] = 10.0;
	}
	else {
		for (int i = 0; i < vec.size(); ++i) vec[i] = 1.0;
	}
}

void GestaltProcessor::NormalizeVec(std::vector< std::vector< float > >& vec) {
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
					vec[i][j] = (vec[i][j] - minValue) / (maxValue - minValue);
				else
					vec[i][j] = 10.0;
	}
	else {
		for (int i = 0; i < vec.size(); ++i)
			for (int j = 0; j < vec[i].size(); ++j)
				vec[i][j] = 1.0;
	}
}

bool GestaltProcessor::CheckMetric(std::vector< std::vector< float > >& vec) {
	for (int i = 0; i < vec.size(); ++i) {
		for (int j = 0; j < vec.size(); ++j) {
			if (vec[i][j] == 0 && i != j) return false;
			if (vec[i][j] != vec[j][i]) return false;
			if (vec[i][j] < 0) return false;
			
			for (int k = 0; k < vec.size(); ++k)
				if (vec[i][k] + vec[k][j] < vec[i][j]) return false;
		}
		if (vec[i][i] != 0) return false;
	}
	return true;
}