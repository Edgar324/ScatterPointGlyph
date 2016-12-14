#ifndef UTILITY_H_
#define UTILITY_H_

#include<vector>
#include <algorithm>
using namespace std;

class CNode;
class ParallelDataset;

class Utility
{
public:
	Utility();
	~Utility();

	// Sort index_two based on index_one
	static void Sort(std::vector<int>& index_one, std::vector<int>& index_two);
    static void Sort(std::vector<float>& index_one, std::vector<int>& index_two);

    static void GridConnection(std::vector<CNode*>& nodes, int w, int h, std::vector<std::vector<bool>>& connecting_status, float& min_edge_length);

	// Delaunay triangulation using vtkDelaunay2D
	// nodes: the input 2D points (Size: N)
	// connecting_status: boolean matrix (N*N)
    static void Triangulation(vector<vector<double>>& pos, vector<vector<bool>>& connecting_status);

	static void GenerateAxisOrder(ParallelDataset* dataset_t, std::vector<int>& axis_order);

    static void GenerateAxisOrder(vector<vector<float>>& values, vector<int>& axis_order);
    static float GetAverageDistance(vector<vector<float>>& pos);
    static float GetDirectAverageDistance(vector<vector<float>>& pos);

    static bool CheckInside(vector<float>& path, float x, float y);
};

#endif