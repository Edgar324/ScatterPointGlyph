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
	void GenerateLayout();
	void GetFinalGlyphPoint(std::vector< int >& point_index);

signals:
	void FinalGlyphUpdated();

private:
	std::vector< std::vector< float > > point_pos_;
	std::vector< float > point_values_;

	/// for kmeans
	int cluster_num_;
	float dis_scale_;
	std::vector< int > cluster_index_;
	std::vector< std::vector< float > > cluster_center_;
	std::vector< int > cluster_node_count_;

	/// for triangulation
	std::vector< std::vector< bool > > connect_status_;

	/// for extracted gestalt
	std::vector< bool > is_node_labeled_;
	std::map< int, int > remaining_node_gestalt_map_;
	std::map< int, int > remaining_gestalt_node_map_;
	std::vector< std::vector< std::vector< int > > > proposal_gestalt_;
	std::vector< float > gestalt_threshold_;

	/// for optimized gestalt
	std::vector< std::vector< int > > optimized_labels_;

	/// for gestalt selection
	std::vector< std::vector< float > > gestalt_info_;

	int final_gestalt_count_;
	std::vector< int > final_label_;

	void KmeansCluster();
	void Triangulation();

	void ProcessGestaltRule(int rule);

	void ExtractProposalGestalt(int rule);
	void ExecLabelOptimization(int rule);

	void Sort(std::vector< std::vector< float > >& value, int begin, int end, int sort_index);
};

#endif