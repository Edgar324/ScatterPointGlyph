#include "perception_solver.h"

PerceptionSolver::PerceptionSolver() {

}

PerceptionSolver::~PerceptionSolver() {

}

void PerceptionSolver::SetPointData(std::vector< std::vector< float > >& pos, std::vector< std::vector< float > >& values) {
	point_pos_ = pos;
	point_values_ = values;

	candidate_points_.resize(point_pos_.size());
	for (size_t i = 0; i < point_pos_.size(); ++i) candidate_points_[i] = i;
}

void PerceptionSolver::ExecRandomSample(int n) {

}

void PerceptionSolver::ExecOptimization() {

}