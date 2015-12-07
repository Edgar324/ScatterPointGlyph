#ifndef CLUSTER_SOLVER_H_
#define CLUSTER_SOLVER_H_

#include <QtCore/QThread>

class ClusterSolver : public QThread
{
	Q_OBJECT

public:
	ClusterSolver();
	~ClusterSolver();

	virtual void GetCluster(int& cluster_count, std::vector< int >& cluster_index);
	void GetFinalClusterPoint(int index, std::vector< int >& point_index);

signals:
	void ClusterUpdated(int index);

protected:
	int final_label_count_;
	std::vector< int > final_label_;

	virtual void run();
};

#endif