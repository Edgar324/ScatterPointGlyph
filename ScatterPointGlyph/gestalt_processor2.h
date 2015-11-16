#ifndef GESTALT_PROCESSOR2_H_
#define GESTALT_PROCESSOR2_H_

#include <vector>
#include <QtCore/QObject>

class ScatterPointDataset;
class GestaltCandidateSet;
class PropertyExtractor;

class GestaltProcessor2 : public QObject
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

	void GenerateCluster(float dis_thresh);
	void GetClusterIndex(int& cluster_count, std::vector< int >& cluster_index);
	void GetFinalGlyphPoint(int index, std::vector< int >& point_index);

signals:
	void FinalGlyphUpdated(int index);

private:
	GestaltCandidateSet* gestalt_candidates_;
	ScatterPointDataset* dataset_;
	std::vector< bool > is_property_on_;
	std::vector< float > property_thresh;
	std::vector< PropertyExtractor* > property_extractors_;

	int extracted_gestalt_count_;
	std::vector< int > final_label_;

	float valid_decreasing_rate_;

	void ExtractLabels();
	void ExtractValidGestalt(int& property_index, int& gestalt_index);

	void SortGestaltInfos(std::vector< GestaltInfo >& infos, int begin, int end, int sort_index);
};

#endif