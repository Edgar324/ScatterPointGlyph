/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "slic_vf.h"
#include <QImage>
#include <QColor>

Slic::Slic()
{

}

Slic::~Slic()
{

}

void Slic::ClearData()
{
    clusters_.clear();
    distances_.clear();
    centers_.clear();
    center_count_.clear();
}

void Slic::InitData()
{
    if (vfd_ == NULL) return;

    for (int i = 0; i < vfd_->width(); i++) {
        vector<int> cr;
        vector<float> dr;
        for (int j = 0; j < vfd_->height(); j++) {
            cr.push_back(-1);
            dr.push_back(DBL_MAX);
        }
        clusters_.push_back(cr);
        distances_.push_back(dr);
    }

    // Initialize the centers and counters.
    for (int i = step_; i < vfd_->width() - step_/2; i += step_) {
        for (int j = step_; j < vfd_->height() - step_/2; j += step_) {
            vector<float> center;
            // Find the local minimum (gradient-wise).
            int yi, xi;
            this->FindLocalMinimum(i, j, xi, yi);
            
            // Generate the center vector.
            center.push_back(xi);
            center.push_back(yi);
            vector<float>& val = vfd_->GetData(xi, yi);

            for (int k = 0; k < val.size(); k++)
                center.push_back(val[k]);
            
            // Append to vector of centers.
            centers_.push_back(center);
            center_count_.push_back(0);
        }
    }
}


float Slic::ComputeDist(int ci, int xi, int yi) {
    vector<float>& val = vfd_->GetData(xi, yi);

    float dc = 0;
    for (int i = 0; i < val.size(); ++i)
        dc += pow(centers_[ci][i + 2] - val[i], 2);

    float ds = sqrt(pow(centers_[ci][0] - xi, 2) + pow(centers_[ci][1] - yi, 2));
    
    return sqrt(pow(dc / nc_, 2) + pow(ds / ns_, 2));
}


void Slic::FindLocalMinimum(int cx, int cy, int& xi, int& yi) {
    float min_grad = DBL_MAX;
    xi = cx;
    yi = cy;
    
    for (int i = cx - 1; i < cx + 2; i++) {
        for (int j = cy - 1; j < cy + 2; j++) {
            vector<float>& c1 = vfd_->GetData(i, j+1);
            vector<float>& c2 = vfd_->GetData(i + 1, j);
            vector<float>& c3 = vfd_->GetData(i, j);

            float i1toi3 = 0;
            for (int k = 0; k < c1.size(); k++) 
                i1toi3 += pow(c1[k] - c3[k], 2);
            i1toi3 = sqrt(i1toi3);

            float i2toi3 = 0;
            for (int k = 0; k < c2.size(); k++)
                i2toi3 += pow(c2[k] - c3[k], 2);
            i2toi3 = sqrt(i2toi3);

            if (i1toi3 + i2toi3 < min_grad) {
                min_grad = i1toi3 + i2toi3;
                xi = i;
                yi = j;
            }
        }
    }
}

void Slic::GenerateSuperPixels(VectorFieldData* vfd, int step, int nc) {
    this->vfd_ = vfd;
    this->step_ = step;
    this->nc_ = nc;
    this->ns_ = step;
    
    /* Clear previous data (if any), and re-initialize it. */
    ClearData();
    InitData();
    
    /* Run EM for 10 iterations (as prescribed by the algorithm). */
    for (int i = 0; i < nr_iterations_; i++) {
        /* Reset distance values. */
        for (int j = 0; j < vfd_->width(); j++) {
            for (int k = 0;k < vfd_->height(); k++) {
                distances_[j][k] = DBL_MAX;
            }
        }

        for (int j = 0; j < (int)centers_.size(); j++) {
            /* Only compare to pixels in a 2 x step by 2 x step region. */
            for (int k = centers_[j][0] - step; k < centers_[j][0] + step; k++) {
                for (int l = centers_[j][1] - step; l < centers_[j][1] + step; l++) {
                
                    if (k >= 0 && k < vfd_->width() && l >= 0 && l < vfd_->height()) {
                        float d = ComputeDist(j, k, l);
                        
                        /* Update cluster allocation if the cluster minimizes the
                           distance. */
                        if (d < distances_[k][l]) {
                            distances_[k][l] = d;
                            clusters_[k][l] = j;
                        }
                    }
                }
            }
        }
        
        /* Clear the center values. */
        for (int j = 0; j < (int) centers_.size(); j++) {
            centers_[j].assign(centers_[j].size(), 0);
            center_count_[j] = 0;
        }
        
        /* Compute the new cluster centers. */
        for (int j = 0; j < vfd_->width(); j++) {
            for (int k = 0; k < vfd_->height(); k++) {
                int c_id = clusters_[j][k];
                
                if (c_id != -1) {
                    vector<float>& val = vfd_->GetData(j, k);
                    
                    centers_[c_id][0] += j;
                    centers_[c_id][1] += k;
                    for (int t = 0; t < val.size(); ++t)
                        centers_[c_id][2 + t] += val[t];
 
                    center_count_[c_id] += 1;
                }
            }
        }

        /* Normalize the clusters. */
        for (int j = 0; j < (int) centers_.size(); j++) {
            for (int k = 0; k < centers_[j].size(); k++)
                centers_[j][k] /= center_count_[j];
        }
    }

    this->SaveContour("contour.jpg");
}

void Slic::CreateConnectivity() {
    int label = 0, adjlabel = 0;
    const int lims = (vfd_->width() * vfd_->height()) / ((int)centers_.size());
    
    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    
    /* Initialize the new cluster matrix. */
    vector<vector<int>> new_clusters;
    for (int i = 0; i < vfd_->width(); i++) { 
        vector<int> nc;
        for (int j = 0; j < vfd_->height(); j++) {
            nc.push_back(-1);
        }
        new_clusters.push_back(nc);
    }

    for (int i = 0; i < vfd_->width(); i++) {
        for (int j = 0; j < vfd_->height(); j++) {
            if (new_clusters[i][j] == -1) {
                vector<Point2D> elements;
                elements.push_back(Point2D(i, j));
            
                /* Find an adjacent label, for possible use later. */
                for (int k = 0; k < 4; k++) {
                    int x = elements[0].x + dx4[k], y = elements[0].y + dy4[k];
                    
                    if (x >= 0 && x < vfd_->width() && y >= 0 && y < vfd_->height()) {
                        if (new_clusters[x][y] >= 0) {
                            adjlabel = new_clusters[x][y];
                        }
                    }
                }
                
                int count = 1;
                for (int c = 0; c < count; c++) {
                    for (int k = 0; k < 4; k++) {
                        int x = elements[c].x + dx4[k], y = elements[c].y + dy4[k];
                        
                        if (x >= 0 && x < vfd_->width() && y >= 0 && y < vfd_->height()) {
                            if (new_clusters[x][y] == -1 && clusters_[i][j] == clusters_[x][y]) {
                                elements.push_back(Point2D(x, y));
                                new_clusters[x][y] = label;
                                count += 1;
                            }
                        }
                    }
                }
                
                /* Use the earlier found adjacent label if a segment size is
                   smaller than a limit. */
                if (count <= lims >> 2) {
                    for (int c = 0; c < count; c++) {
                        new_clusters[elements[c].x][elements[c].y] = adjlabel;
                    }
                    label -= 1;
                }
                label += 1;
            }
        }
    }
}

void Slic::SaveContour(const char* file_name) {
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	
	/* Initialize the contour vector and the matrix detailing whether a pixel
	 * is already taken to be a contour. */
	vector<Point2D> contours;
	vector<vector<bool>> istaken;
	for (int i = 0; i < vfd_->width(); i++) { 
        vector<bool> nb;
        for (int j = 0; j < vfd_->height(); j++) {
            nb.push_back(false);
        }
        istaken.push_back(nb);
    }
    
    /* Go through all the pixels. */
    for (int i = 0; i < vfd_->width(); i++) {
        for (int j = 0; j < vfd_->height(); j++) {
            int nr_p = 0;
            
            /* Compare the pixel to its 8 neighbours. */
            for (int k = 0; k < 8; k++) {
                int x = i + dx8[k], y = j + dy8[k];
                
                if (x >= 0 && x < vfd_->width() && y >= 0 && y < vfd_->height()) {
                    if (istaken[x][y] == false && clusters_[i][j] != clusters_[x][y]) {
                        nr_p += 1;
                    }
                }
            }
            
            /* Add the pixel to the contour list if desired. */
            if (nr_p >= 2) {
                contours.push_back(Point2D(i,j));
                istaken[i][j] = true;
            }
        }
    }
    
    QImage image(vfd_->width(), vfd_->height(), QImage::Format_RGB32);
    for (int i = 0; i < vfd_->height(); ++i)
        for (int j = 0; j < vfd_->width(); ++j) {
            QColor color;
            vector<float>& val = vfd_->GetData(j, i);
            color.setRedF(val[2]);
            color.setGreenF(val[2]);
            color.setBlueF(val[2]);
            image.setPixel(j, i, color.rgb());
        }

    /* Draw the contour pixels. */
    QColor contour_color;
    contour_color.setRgb(255, 0, 0);
    for (int i = 0; i < (int)contours.size(); i++) {
        image.setPixel(contours[i].x, contours[i].y, contour_color.rgb());
    }

    image.save(file_name);
}