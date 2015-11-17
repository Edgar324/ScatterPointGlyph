#ifndef GESTALT_PROCESSOR2_H_
#define GESTALT_PROCESSOR2_H_

#include <vector>
#include "cluster_solver.h"

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

	struct GestaltInfo {
		//energy, decreasing_rate;
		float values[2];
		int property_index;
		int gestalt_index;
	};

	void SetData(ScatterPointDataset* data);
	void SetPropertyOn(GestaltProperty property);
	void SetPropertyOff(GestaltProperty property);
	void SetThreshold(GestaltProperty peroperty, float thresh);
	void SetDisThreshold(float dis_thresh);

signals:

protected:
	virtual void run();

private:
	GestaltCandidateSet* gestalt_candidates_;
	ScatterPointDataset* dataset_;
	std::vector< bool > is_property_on_;
	std::vector< float > property_thresh_;
	std::vector< PropertyExtractor* > property_extractors_;

	float valid_decreasing_rate_;
	float dis_threshold_;
	int labeled_site_count_;

	void GenerateCluster(float dis_thresh);

	void ExtractLabels();
	void ExtractValidGestalt(int& property_index, int& gestalt_index);

	void SortGestaltInfos(std::vector< GestaltInfo >& infos, int begin, int end, int sort_index);
};

#endif