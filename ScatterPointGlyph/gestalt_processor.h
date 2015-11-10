#ifndef GESTALT_PROCESSOR_H_
#define GESTALT_PROCESSOR_H_

#include <vector>
#include <triangle.h>

class GestaltProcessor
{
public:
	GestaltProcessor();
	~GestaltProcessor();

	void SetData(std::vector< std::vector< float > >& pos, std::vector< std::vector< float > >& values, std::vector< float >& weight);

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
	struct triangulateio triangle_in, triangle_mid, triangle_out, triangle_vorout;
	std::vector< std::vector< bool > > connect_status_;

	/// for extracted gestalt
	std::vector< bool > is_node_labeled_;
	std::vector< std::vector< std::vector< int > > > proposal_gestalt_;
	std::vector< float > gestalt_threshold_;

	/// for optimized gestalt
	std::vector< std::vector< int > > optimized_labels_;

	/// for gestalt selection
	std::vector< std::vector< float > > gestalt_energy_;
	int final_gestalt_count_;
	std::vector< int > final_label_;

	void KmeansCluster();
	void Triangulation();
	void GenerateLayout();

	void ProcessGestaltRule(int rule);

	void ExtractProposalGestalt(int rule);
	void ExecLabelOptimization(int rule);
};

#endif