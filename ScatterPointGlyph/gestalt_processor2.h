#ifndef GESTALT_PROCESSOR2_H_
#define GESTALT_PROCESSOR2_H_

#include <vector>
#include "cluster_solver.h"
#include "tree_common.h"

class ScatterPointDataset;
class GestaltCandidateSet;
class PropertyExtractor;

class GestaltProcessor2 : public ClusterSolver
{
	Q_OBJECT

public:
	GestaltProcessor2();
	~GestaltProcessor2();

	enum GestaltProperty {
		PROXIMITY = 0x0,
		SIMILARITY,
		CONTINUITY,
		COMMON_FATE,
		CLOSURE
	};

	void SetData(std::vector< CLeaf* >& nodes, std::vector< CNode* >& clusters,
		std::vector< int >& node_index, std::vector< std::vector< bool > >& node_connecting_status,
		std::vector< float >& var_weights);

	void SetPropertyOn(GestaltProperty property);
	void SetPropertyOff(GestaltProperty property);
	void SetThreshold(GestaltProperty property, float thresh);
	void SetDisThreshold(float dis_thresh);

	void GetResultLabel(std::vector< int >& result);

	virtual void GenerateCluster();
	virtual void GetCluster(int& cluster_count, std::vector< int >& cluster_index);
	virtual void GetCluster(float fitness, std::vector< int >& cluster_index);

protected:
	virtual void run();

private:
	GestaltCandidateSet* gestalt_candidates_;

	float dis_threshold_;

	std::vector< bool > is_property_on_;
	std::vector< float > property_thresh_;
	std::vector< float > property_weight_;
	std::vector< PropertyExtractor* > property_extractors_;

	std::vector< int > result_label_;

	void ExtractLabels();
};

#endif