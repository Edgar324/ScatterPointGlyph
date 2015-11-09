#ifndef HIER_SOVLER_H_
#define HIER_SOLVER_H_

#include <QtCore/QThread>
#include <vector>

class HierSolver : public QThread
{
	Q_OBJECT

public:
	HierSolver();
	~HierSolver();

	void SetData(std::vector< std::vector< float > >& value, std::vector< float >& weight);
	void SetInitialClusterNum(int num);
	void SetClusterNum(int num);
	int GetCurrentClusterNum() { return current_cluster_num_; }
	void GetPreClusterIndex(int& cluster_one, int& cluster_two);

	void GetGlyphData(std::vector< float >& glyph_pos, std::vector< std::vector< float > >& glyph_values);

signals:
	void ClusterUpdated();
	void CombinedClusterChanged(int, int);

protected:
	void run();

private:
	std::vector< std::vector< float > > value_;
	std::vector< float > weight_;
	std::vector< int > cluster_index_;
	std::vector< std::vector< float > > cluster_center_;
	std::vector< int > cluster_node_count_;

	int cluster_num_;

	int current_cluster_num_;
	int pre_cluster_index_[2];

	void KmeansCluster();
};

#endif