#ifndef PARALLEL_COORDINATE_H_
#define PARALLEL_COORDINATE_H_

#include "gl/glew.h"
#include <QtOpenGL/QGLWidget>
#include <vector>
#include <iostream>
#include "parallel_dataset.h"

class ParallelCoordinate : public QGLWidget {
    Q_OBJECT

public:
    ParallelCoordinate();
    ~ParallelCoordinate();

    enum PcpPlotStatus{
        PLOT_OK				= 0x000001,
        GL_INIT_ERROR	    = 0x000010
    };

    void SetData(ParallelDataset* dataset_t);
	void SetAxisOrder(std::vector<int>& axis_order);
	void SetHighlightAxis(int var_index);

public slots:
    void OnDataChanged();

signals:
	void HighlightVarChanged(int);

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
	void mouseDoubleClickEvent(QMouseEvent *);
	void mousePressEvent(QMouseEvent *);

private:
    ParallelDataset* dataset_;

    float x_border_, y_border_, axis_name_height_, range_text_height_, axis_width_, weight_circle_radius_;
    float subset_rect_width_, subset_rect_height_, subset_rect_text_width_, subset_rect_y_value_;

    float axis_top_y_value_, axis_bottom_y_value_, axis_y_size_;
    std::vector<float> axis_x_pos_values_;
    float axis_name_y_value_, range_text_top_y_value_, range_text_bottom_y_value_;

    float weight_circle_center_y_value_;

    GLuint setting_texture_;
    float icon_width_, icon_height_;

	int highlight_var_index_;

    void UpdateViewLayoutParameters();

    void PaintSettingIcon();
    void PaintCoordinate();
    void PaintLines();
    void PaintText();
    void PaintSubsetIdentifyItems();
    void PaintWeightCircles();
	void PaintGaussianCurve();

    void UpdateView();
};

#endif