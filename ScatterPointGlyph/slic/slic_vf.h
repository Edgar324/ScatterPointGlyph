/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include <vector>
using namespace std;

struct Point2D {
    int x, y;

    Point2D(int v1, int v2) : x(v1), y(v2) {}
};

class VectorFieldData
{
public:
    VectorFieldData(vector<vector<vector<double>>>& v)
        : values_(v) 
    {
        if (values_.size() > 0 && values_[0].size() > 0 && values_[0][0].size() > 0) {
            h_ = values_.size();
            w_ = values_[0].size();
            v_num_ = values_[0][0].size();
        }
    }

    ~VectorFieldData() {}

    int width() { return w_; }
    int height() { return h_; }
    int v_num() { return v_num_; }
    vector<double>& GetData(int x, int y) { return values_[y][x]; }

private:
    int w_ = 0, h_ = 0, v_num_ = 0;
    vector<vector<vector<double>>> values_;
};

class Slic
{
public:
    Slic();
    ~Slic();

    void GenerateSuperPixels(VectorFieldData* vfd, int step, int nc);
    vector<vector<int>>& GetClusters() { return clusters_; }
    vector<vector<double>>& GetCenters() { return centers_; }

private:
    VectorFieldData* vfd_ = NULL;

    const int nr_iterations_ = 10;

    // The cluster assignments and distance values for each pixel.
    vector<vector<int>> clusters_;
    vector<vector<double>> distances_;

    // The xy and vector values of the centers.
    vector<vector<double>> centers_;
    // The number of occurences of each center.
    vector<int> center_count_;

    // The step size per cluster, and the colour (nc) and distance (ns) parameters.
    int step_, nc_, ns_;

    // Compute the distance between a center and an individual pixel
    double ComputeDist(int ci, int xi, int yi);
    
    // Find the pixel with the lowest gradient in a 3x3 surrounding.
    void FindLocalMinimum(int cx, int cy, int& xi, int& yi);

    void CreateConnectivity();
        
    // Remove and initialize the 2d vectors.
    void ClearData();
    void InitData();

    void SaveContour(const char* file_name);
};