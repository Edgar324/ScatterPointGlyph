#ifndef UTILITY_H_
#define UTILITY_H_

#include<vector>

class CNode;

class Utility
{
public:
	Utility();
	~Utility();

	// Sort index_two based on index_one
	static void Sort(std::vector<int>& index_one, std::vector<int>& index_two);

	// Delaunay triangulation using vtkDelaunay2D
	// nodes: the input 2D points (Size: N)
	// connecting_status: boolean matrix (N*N)
	// min_edge_length: the minimum distance between any two points
	static void VtkTriangulation(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status, float& min_edge_length);
};

#endif