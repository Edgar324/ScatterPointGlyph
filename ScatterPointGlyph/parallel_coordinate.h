#ifndef PARALLEL_COORDINATE_H_
#define PARALLEL_COORDINATE_H_

#include "gl/glew.h"
#include <QtGui/QColor>
#include <QtOpenGL/QGLWidget>
#include <vector>
#include <iostream>

class ParallelRecord {
public:
    ParallelRecord() { }        
    ~ParallelRecord() { }

    std::vector< float > values;
};

class ParallelDataset {
public:
    ParallelDataset();
    ~ParallelDataset();

    bool CompleteInput();
	void UpdateGaussian();
    bool ClearData();

    // attributes which must be set
    std::vector< QString > subset_names;
    std::vector< std::vector< ParallelRecord* > > subset_records;
    std::vector< std::vector< QString > > axis_anchors;
    std::vector< QString > axis_names;
	std::vector< QColor > subset_colors;


    // attributes which can be set automatically
	std::vector< std::vector< float > > var_centers;
	std::vector< std::vector< float > > var_width;
    std::vector< std::vector< QColor > > record_color;
    std::vector< std::vector< bool > > is_record_selected;

    std::vector< bool > is_subset_visible;
    std::vector< float > subset_opacity;
    std::vector< bool > is_axis_selected;
    std::vector< int > mapped_axis;

	bool is_axis_weight_enabled;
	std::vector< float > axis_weights;

    bool is_cluster_enabled;
    std::vector< std::vector< ParallelRecord* > > cluster_centers;

    bool is_range_filter_enabled;
    std::vector< std::vector< float > > axis_value_filter_range;

    bool is_edge_bundling_enabled;
    bool is_correlation_analysis_enabled;
	bool is_gaussian_enabled;

	bool is_updating;
};

class ParallelCoordinate : public QGLWidget {
    Q_OBJECT

public:
    ParallelCoordinate();
    ~ParallelCoordinate();

    enum PcpPlotStatus{
        PLOT_OK				= 0x000001,
        GL_INIT_ERROR	    = 0x000010
    };


    void SetDataset(ParallelDataset* dataset_t);

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

private:
    ParallelDataset* dataset_;

    float x_border_, y_border_, axis_name_height_, range_text_height_, axis_width_, weight_circle_radius_;
    float subset_rect_width_, subset_rect_height_, subset_rect_text_width_, subset_rect_y_value_;

    float axis_top_y_value_, axis_bottom_y_value_, axis_y_size_;
    std::vector< float > axis_x_pos_values_;
    float axis_name_y_value_, range_text_top_y_value_, range_text_bottom_y_value_;

    float weight_circle_center_y_value_;

    GLuint setting_texture_;
    float icon_width_, icon_height_;

    void UpdateViewLayoutParameters();

    void PaintSettingIcon();
    void PaintCoordinate();
    void PaintLines();
    void PaintText();
    void PaintSubsetIdentifyItems();
    void PaintWeightCircles();

	void PaintGaussianCurve();
};

#endif