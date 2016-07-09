#ifndef VOLUME_RENDER_WIDGET_H_
#define VOLUME_RENDER_WIDGET_H_

#include "GL/glew.h"
#include "glm/glm.hpp"
#include <QtOpenGL/QGLWidget>
#include <QtCore/QTimer>
#include <vector>

class ArcBall;
class GlProgram;

class VolumeRenderWidget : public QGLWidget
{
    Q_OBJECT

public:
    VolumeRenderWidget();
    ~VolumeRenderWidget();

    void SetData(int* sizes_t, float* spacings_t, GLenum data_format_t, void* data_t);
    void SetBackgroundData(int w, int h, std::vector< float > image_data);
    void SetTransferFunction(float* values, int entryNumber);
    void SetViewWindow(float win_center_t, float win_width_t);
    void SetSampleStep(float sample_step_t);

    void OptimizeResult();

    float min_data_value() { return min_data_value_; }
    float max_data_value() { return max_data_value_; }

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

private:
    int sizes_[3];
    float spacings_[3];
    GLenum data_format_;
    void* data_;
    float min_data_value_, max_data_value_;

    std::vector< float > tf_values_;
    std::vector< float > random_texture_values_;
    std::vector< float > background_image_data_;

    float zoom_factor_;
    ArcBall* rotation_ball_;

    GLuint VAO_, VBO_, EBO_;
    GLuint FBO_;

    GlProgram* back_face_prog_;
    GlProgram* render_volume_prog_;
    GlProgram* render_background_prog_;

    GLuint transfer_function_tex_;
    GLuint back_face_tex_;
    GLuint random_tex_;
    GLuint volume_tex_;
    GLuint depth_render_buffer_;
    GLuint optimization_result_tex_[2];
    GLuint background_tex_;

    float sample_step_;
    float volume_bias_, volume_scale_;

    glm::vec3 back_ground_color_, light_dir_, light_half_vec_;
    glm::vec4 light_diffuse_, light_ambient_, light_specular_;
    float shininess_;

    int optimize_run_;
    bool is_optimizing_;
    std::vector< float > optimization_pixel_values_, rendering_result_;

    QTimer* optimizing_timer_;

    float map_left_, map_right_, map_bottom_, map_top_;

    void InitShaderProgram();
    void InitArrayBuffer();
    void InitFrameBuffer();
    void InitTextures();

    void UpdateVolumeTexture();
    void UpdateTransferTexture();
    void UpdateBackFaceTexture();
    void UpdateRandomTexture();
    void UpdateOptimizationResultTexture();

    void UpdateBackFacePara();
    void UpdateRenderVolumePara();
    void UpdateBackgroundPara();

    template < typename T >
    void GetMinMaxValue(unsigned int data_size, T* data, float& min_value, float& max_value);

    void OutputErrorLogAndStopRendering(const char* message);
    void RestartOptimizingTimer();

    private slots:
        void OnOptimizingTimerTimeOut();
};

#endif