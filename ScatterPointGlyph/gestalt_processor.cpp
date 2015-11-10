#include "gestalt_processor.h"
#include <time.h>
#include <queue>
#include <iostream>
#include "./gco-v3.0/GCoptimization.h"

GestaltProcessor::GestaltProcessor() 
	: cluster_num_(200), dis_scale_(0.5){

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

void GestaltProcessor::GenerateLayout() {
	// step one: kmeans to generate a series of candidate sites
	KmeansCluster();

	// step two: generate connecting graph using delaunay triangulation
	Triangulation();

	is_node_labeled_.resize(cluster_num_);
	is_node_labeled_.assign(cluster_num_, false);

	final_gestalt_count_ = -1;
	final_label_.resize(cluster_num_);
	final_label_.assign(cluster_num_, -1);

	bool is_all_labeled = false;
	while (!is_all_labeled) {
		// step three: generate gestalts based on different rules
		// 1) Similarity   2) Proximity  3) Continuity  4) closure
		// step four: generate optimization result based on the nodes
		ProcessGestaltRule(0);
		ProcessGestaltRule(1);

		// step five: remove the strongest gestalt, goto step four until all nodes are removed
		int rule_index = -1, gestalt_index = -1;
		int max_energy = -1e10;
		for (int i = 0; i < gestalt_energy_.size(); ++i)
			for (int j = 0; j < gestalt_energy_[i].size(); ++j)
				if (gestalt_energy_[i][j] > max_energy) {
					max_energy = gestalt_energy_[i][j];
					rule_index = i;
					gestalt_index = j;
				}

		final_gestalt_count_++;
		for (int i = 0; i < proposal_gestalt_[rule_index][gestalt_index].size(); ++i) {
			is_node_labeled_[proposal_gestalt_[rule_index][gestalt_index][i]] = true;
			final_label_[proposal_gestalt_[rule_index][gestalt_index][i]] = final_gestalt_count_;
		}
		is_all_labeled = true;
		for (int i = 0; i < is_node_labeled_.size(); ++i)
			if (!is_node_labeled_[i]) {
				is_all_labeled = false;
				break;
			}
	}
}

void GestaltProcessor::KmeansCluster() {
	/// Generate random center
	cluster_center_.resize(cluster_num_);
	cluster_node_count_.resize(cluster_num_);

	std::vector< bool > is_selected;
	is_selected.resize(point_values_.size(), false);
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
}

void GestaltProcessor::Triangulation() {
	triangle_in.numberofpoints = cluster_num_;
	triangle_in.numberofpointattributes = 0;
	triangle_in.pointlist = (float *)malloc(triangle_in.numberofpoints * 2 * sizeof(float));
	for (int i = 0; i < cluster_num_; ++i){
		triangle_in.pointlist[i * 2] = cluster_center_[i][0];
		triangle_in.pointlist[i * 2 + 1] = cluster_center_[i][1];
	}
	triangle_in.pointattributelist = NULL;
	triangle_in.pointmarkerlist = NULL;

	triangle_in.numberofsegments = 0;
	triangle_in.numberofholes = 0;
	triangle_in.numberofregions = 1;
	triangle_in.regionlist = (float *)malloc(triangle_in.numberofregions * 4 * sizeof(float));
	triangle_in.regionlist[0] = 0.5;
	triangle_in.regionlist[1] = 5.0;
	triangle_in.regionlist[2] = 7.0;            /* Regional attribute (for whole mesh). */
	triangle_in.regionlist[3] = 0.1;          /* Area constraint that will not be used. */

	triangle_mid.pointlist = (float *)NULL;            /* Not needed if -N switch used. */
	/* Not needed if -N switch used or number of point attributes is zero: */
	triangle_mid.pointattributelist = (float *)NULL;
	triangle_mid.pointmarkerlist = (int *)NULL; /* Not needed if -N or -B switch used. */
	triangle_mid.trianglelist = (int *)NULL;          /* Not needed if -E switch used. */
	/* Not needed if -E switch used or number of triangle attributes is zero: */
	triangle_mid.triangleattributelist = (float *)NULL;
	triangle_mid.neighborlist = (int *)NULL;         /* Needed only if -n switch used. */
	/* Needed only if segments are output (-p or -c) and -P not used: */
	triangle_mid.segmentlist = (int *)NULL;
	/* Needed only if segments are output (-p or -c) and -P and -B not used: */
	triangle_mid.segmentmarkerlist = (int *)NULL;
	triangle_mid.edgelist = (int *)NULL;             /* Needed only if -e switch used. */
	triangle_mid.edgemarkerlist = (int *)NULL;   /* Needed if -e used and -B not used. */

	triangle_vorout.pointlist = (float *)NULL;        /* Needed only if -v switch used. */
	/* Needed only if -v switch used and number of attributes is not zero: */
	triangle_vorout.pointattributelist = (float *)NULL;
	triangle_vorout.edgelist = (int *)NULL;          /* Needed only if -v switch used. */
	triangle_vorout.normlist = (float *)NULL;         /* Needed only if -v switch used. */

	triangulate("pczAevn", &triangle_in, &triangle_mid, &triangle_vorout);

	triangle_out.pointlist = (float *)NULL;            /* Not needed if -N switch used. */
	/* Not needed if -N switch used or number of attributes is zero: */
	triangle_out.pointattributelist = (float *)NULL;
	triangle_out.trianglelist = (int *)NULL;          /* Not needed if -E switch used. */
	/* Not needed if -E switch used or number of triangle attributes is zero: */
	triangle_out.triangleattributelist = (float *)NULL;

	/* Refine the triangulation according to the attached */
	/*   triangle area constraints.                       */

	triangulate("pzBP", &triangle_mid, &triangle_out, (struct triangulateio *) NULL);

	connect_status_.resize(cluster_num_);
	for (int i = 0; i < connect_status_.size(); ++i){
		connect_status_[i].resize(cluster_num_, false);
	}
	for (int i = 0; i < triangle_out.numberoftriangles; ++i){
		int id1 = triangle_out.trianglelist[i * 3];
		int id2 = triangle_out.trianglelist[i * 3 + 1];
		int id3 = triangle_out.trianglelist[i * 3 + 2];
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
	ExecLabelOptimization(rule);

	// step three: generate gestalt energy
}

void GestaltProcessor::ExtractProposalGestalt(int rule){
	proposal_gestalt_[rule].clear();

	std::vector< float > distance;
	distance.resize(cluster_num_);
	std::vector< bool > is_reached;
	is_reached.resize(cluster_num_);

	for (int i = 0; i < cluster_num_; ++i) {
		distance.assign(distance.size(), 1e10);
		is_reached.assign(cluster_num_, false);
		distance[i] = 0;

		for (int j = 0; j < cluster_num_ - 1; ++j) {
			float min_dist = 1e11;
			int min_index = -1;
			for (int k = 0; k < cluster_num_; ++k)
				if (!is_reached[k] && distance[k] < min_dist) {
					min_dist = distance[k];
					min_index = k;
				}
			is_reached[min_index] = true;

			for (int k = 0; k < cluster_num_; ++k)  
				if (!is_reached[k] && connect_status_[min_index][k]) {
					float temp_dis = 0;
					switch (rule)
					{
					case 0:
						temp_dis += abs(cluster_center_[j][2] - cluster_center_[k][2]);							
						break;
					case 1:
						temp_dis += sqrt(pow(cluster_center_[j][0] - cluster_center_[k][0], 2) + pow(cluster_center_[j][1] - cluster_center_[k][1], 2));
						break;
					case 2:
						break;
					case 3:
						break;
					default:
						break;
					}

					if (distance[k] > distance[min_index] + temp_dis) distance[k] = distance[min_index] + temp_dis;
				}
		}

		std::vector< int > gestalt;
		for (int i = 0; i < distance.size(); ++i)
			if (distance[i] < gestalt_threshold_[rule]) gestalt.push_back(i);

		/// TODO: filter out identical gestalt
		proposal_gestalt_[rule].push_back(gestalt);
	}
}

void GestaltProcessor::ExecLabelOptimization(int rule) {
	try {
		// Step 1: construct class
		GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(cluster_num_, cluster_num_);

		// Step 2: construct data cost, smooth cost and label cost for each label
		std::vector< std::vector< float > > data_cost, smooth_cost;
		std::vector< float > label_cost;
		data_cost.resize(cluster_num_);
		smooth_cost.resize(cluster_num_);
		label_cost.resize(cluster_num_, -1);
		for (int i = 0; i < cluster_num_; ++i) {
			data_cost[i].resize(cluster_num_, -1);
			smooth_cost[i].resize(cluster_num_, -1);
		}

		// label cost
		std::vector< std::vector< float > > gestalt_average;
		std::vector< int > gestalt_node_count;
		std::vector< float > gestalt_variance;
		gestalt_average.resize(cluster_num_);
		gestalt_node_count.resize(cluster_num_, 0);
		gestalt_variance.resize(cluster_num_, 0);
		gestalt_average.resize(cluster_num_);
		int average_length = 1;
		if (rule == 1) average_length = 2;
		for (int i = 0; i < cluster_num_; ++i) gestalt_average[i].resize(average_length, 0);

		for (int i = 0; i < cluster_num_; ++i) {
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
		for (int i = 0; i < cluster_num_; ++i)
			for (int j = 0; j < average_length; ++j)
				gestalt_average[i][j] /= gestalt_node_count[i];

		for (int i = 0; i < cluster_num_; ++i) {
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
		for (int i = 0; i < cluster_num_; ++i)
			gestalt_variance[i] = sqrt(gestalt_variance[i] / gestalt_node_count[i]);
		for (int i = 0; i < label_cost.size(); ++i) label_cost[i] = gestalt_variance[i];

		// data cost
		for (int i = 0; i < proposal_gestalt_[rule].size(); ++i) {
			for (int j = 0; j < proposal_gestalt_[rule][i].size(); ++j) {
				int temp_index = proposal_gestalt_[rule][i][j];
				switch (rule){
				case 0:
					data_cost[temp_index][i] = abs(cluster_center_[temp_index][2] - gestalt_average[i][0]);
					break;
				case 1:
					data_cost[temp_index][i] = sqrt(pow(cluster_center_[temp_index][0] - gestalt_average[i][0], 2) + pow(cluster_center_[temp_index][1] - gestalt_average[i][1], 2));
					break;
				default:
					break;
				}
			}
		}
		for (int i = 0; i < cluster_num_; ++i){
			data_cost[i][i] = 0;
			for (int j = 0; j < cluster_num_; ++j)
				if (data_cost[i][j] < 0) data_cost[i][j] = 1e10;
		}

		// smooth cost
		for (int i = 0; i < cluster_num_; ++i)
			for (int j = 0; j < cluster_num_; ++j)
				if (connect_status_[i][j]) {
					smooth_cost[i][j] = 1.0;
					smooth_cost[j][i] = smooth_cost[i][j];
				}
				else
					smooth_cost[i][j] = 0.5;


		// Step 3: multi-label graph cut
		for (int i = 0; i < cluster_num_ - 1; ++i)
			for (int j = i + 1; j < cluster_num_; ++j)
				if (connect_status_[i][j]) gc->setNeighbors(i, j);
		for (int i = 0; i < cluster_num_; ++i) {
			for (int j = 0; j < cluster_num_; ++j) {
				gc->setDataCost(i, j, data_cost[i][j]);
				gc->setSmoothCost(i, j, smooth_cost[i][j]);
			}
			gc->setLabelCost(label_cost[i]);
		}

		std::cout << "\nBefore optimization energy is " << gc->compute_energy() << std::endl;
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		std::cout << "After optimization energy is " << gc->compute_energy() << std::endl;

		optimized_labels_[rule].resize(cluster_num_);
		for (int i = 0; i < cluster_num_; i++)
			optimized_labels_[rule][i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e) {
		e.Report();
	}
}