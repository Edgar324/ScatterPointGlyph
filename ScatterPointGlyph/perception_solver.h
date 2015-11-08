#ifndef PERCEPTION_SOLVER_H_
#define PERCEPTION_SOLVER_H_

#include <vector>

class PerceptionSolver
{
public:
	PerceptionSolver();
	~PerceptionSolver();

	void SetPointData(std::vector< std::vector< float > >& pos, std::vector< std::vector< float > >& values);

	void ExecRandomSample(int n);

	void ExecOptimization();

private:
	std::vector< std::vector< float > > point_pos_;
	std::vector< std::vector< float > > point_values_;

	std::vector< int > candidate_points_;
};

#endif