#ifndef GESTALT_PROCESSOR_H_
#define GESTALT_PROCESSOR_H_

#include <vector>
#include <QtCore/QObject>

class GestaltProcessor : public QObject
{
	Q_OBJECT

public:
	GestaltProcessor();
	~GestaltProcessor();

	void SetData(std::vector< std::vector< float > >& pos, std::vector< std::vector< float > >& values, std::vector< float >& weight);
	void SetGestaltThreshold(std::vector< float >& threshold);
	void SetDistanceThreshold(float thresh);
	void GenerateLayout();
	void GetFinalGlyphPoint(std::vector< int >& point_index);
	void GetKmeansClusterResult(std::vector< std::vector< int > >& index);

signals:
	void FinalGlyphUpdated();
	void KMeansClusterFinished();

private:
	std::vector< std::vector< float > > point_pos_;
	std::vector< float > point_values_;

	/// for kmeans
	bool is_sample_needed_;
	int cluster_num_;
	float dis_scale_;
	std::vector< int > cluster_index_;
	std::vector< std::vector< float > > cluster_center_;
	std::vector< int > cluster_node_count_;

	/// for triangulation
	std::vector< std::vector< bool > > connect_status_;

	/// for extracted gestalt
	std::vector< bool > is_node_labeled_;
	std::vector< std::vector< std::vector< int > > > proposal_gestalt_;
	std::vector< float > gestalt_threshold_;
	float distance_threshold_;
	std::vector< std::vector< float > > node_distance_;

	/// for optimized gestalt
	std::vector< std::vector< int > > optimized_labels_;

	/// for gestalt selection
	std::vector< std::vector< float > > gestalt_info_;

	int final_gestalt_count_;
	std::vector< int > final_label_;

	void KmeansCluster();
	void Triangulation();

	void ProcessGestaltRule(int rule);

	void ExtractNodeDistance();
	void ExtractProposalGestalt(int rule);
	void ExecLabelOptimization(int rule);
	void UpdateDecreasingRate(int rule);

	void Sort(std::vector< std::vector< float > >& value, int begin, int end, int sort_index);
	void NormalizeVec(std::vector< float >& vec);
	void NormalizeVec(std::vector< std::vector< float > >& vec);
	bool CheckMetric(std::vector< std::vector< float > >& vec);
};

#endif