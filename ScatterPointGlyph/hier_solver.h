#ifndef HIER_SOVLER_H_
#define HIER_SOLVER_H_

#include "cluster_solver.h"
#include <vector>

class ScatterPointDataset;

class HierSolver : public ClusterSolver
{
	Q_OBJECT

public:
	HierSolver();
	~HierSolver();

	virtual void GetClusterIndex(int& cluster_count, std::vector< int >& cluster_index) {}
	void SetData(ScatterPointDataset* data);
	void SetExpectedClusterNum(int num);

protected:
	virtual void run();

private:
	ScatterPointDataset* dataset_;

	int expected_cluster_num_;
	int current_cluster_num_;
	std::vector< std::vector< float > > cluster_center_;
	std::vector< int > cluster_node_count_;
	std::vector< bool > is_label_used_;

	std::vector< std::vector< bool > > cluster_connecting_status_;

	float dis_scale_;

	void Triangulation();
};

#endif