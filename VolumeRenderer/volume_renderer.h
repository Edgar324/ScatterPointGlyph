#ifndef VOLUME_RENDERER_H
#define VOLUME_RENDERER_H

#include "GL/glew.h"
#include "volume_renderer_global.h"
#include <vector>
using namespace std;
#include <QWidget>
#include <QColor>

class TransferFunction1DWidget;
class VolumeRenderWidget;
class QSlider;
class QImage;

class VOLUMERENDERER_EXPORT VolumeRenderer : public QWidget
{
    Q_OBJECT

public:
    VolumeRenderer();
    ~VolumeRenderer();

    void SetData(int* sizes_t, float* spacings_t, GLenum data_format_t, void* data_t);
    void SetData(int* sizes_t, float* spacings_t, unsigned char* data_t, vector<QColor>& colors);
    QImage* GetRenderingImage();

private:
    TransferFunction1DWidget* transfer_function_1d_widget_;
    VolumeRenderWidget* volume_render_widget_;
    QSlider* win_center_slider;
    QSlider* win_width_slider;

    const int TF_ENTRY_NUM = 256;

    void InitializeWidget();
        
    private slots:
        void OnTransferFunctionChanged();
        void OnWinCenterChanged(int);
        void OnWinWidthChanged(int);
        void OnOptimizationTriggered();
};

#endif // VOLUME_RENDERER_H
