#include "cluster_solver.h"

ClusterSolver::ClusterSolver() 
	: final_label_count_(0) {

}

ClusterSolver::~ClusterSolver() {

}

void ClusterSolver::GetClusterIndex(int& cluster_count, std::vector< int >& cluster_index) {
	cluster_count = final_label_count_;
	cluster_index = final_label_;
}

void ClusterSolver::GetFinalClusterPoint(int index, std::vector< int >& point_index) {
	point_index.clear();
	for (int i = 0; i < final_label_.size(); ++i)
		if (final_label_[i] == index) point_index.push_back(i);
}

void ClusterSolver::run() {

}