#include "volume_render_widget.h"
#include <iostream>
#include <limits.h>
using namespace std;
#include <QMessageBox>
#include <QMouseEvent>
#include <QWheelEvent>
#include "glshader.h"
#include "glprogram.h"
#include "arcball.h"

VolumeRenderWidget::VolumeRenderWidget()
    : data_(0), zoom_factor_(1.0),
      VAO_(-1), VBO_(-1), EBO_(-1), FBO_(-1), back_face_prog_(0), render_volume_prog_(0), sample_step_(0.005), 
      transfer_function_tex_(-1), back_face_tex_(-1), volume_tex_(-1), depth_render_buffer_(-1), is_optimizing_(false){

    this->setMinimumHeight(300);
    this->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    map_left_ = 70;
    map_right_ = 138;
    map_bottom_ = 15;
    map_top_ = 60;

    rotation_ball_ = new ArcBall;
    back_ground_color_[0] = 1.0;
    back_ground_color_[1] = 1.0;
    back_ground_color_[2] = 1.0;

    light_dir_[0] = 1;
    light_dir_[1] = 1;
    light_dir_[2] = 1;
    light_dir_ = glm::normalize(light_dir_);

    light_ambient_[0] = 1.0;
    light_ambient_[1] = 1.0;
    light_ambient_[2] = 1.0;
    light_ambient_[3] = 1.0;

    light_diffuse_[0] = 0.1;
    light_diffuse_[1] = 0.1;
    light_diffuse_[2] = 0.1;
    light_diffuse_[3] = 0.1;

    light_specular_[0] = 1;
    light_specular_[1] = 1;
    light_specular_[2] = 1;
    light_specular_[3] = 1;

    shininess_ = 5;

    optimize_run_ = 0;
    optimizing_timer_ = new QTimer;
    optimizing_timer_->setSingleShot(true);
    connect(optimizing_timer_, SIGNAL(timeout()), this, SLOT(OnOptimizingTimerTimeOut()));
}

VolumeRenderWidget::~VolumeRenderWidget(){

}

void VolumeRenderWidget::SetData(int* sizes_t, float* spacings_t, GLenum data_format_t, void* data_t){
    int data_size = sizes_t[0] * sizes_t[1] * sizes_t[2];
    switch (data_format_t){
    case GL_UNSIGNED_BYTE:
        GetMinMaxValue(data_size, (unsigned char*)data_t, min_data_value_, max_data_value_);
        break;
    case GL_SHORT:
        GetMinMaxValue(data_size, (short*)data_t, min_data_value_, max_data_value_);
        break;
    case GL_FLOAT:{
        GetMinMaxValue(data_size, (float*)data_t, min_data_value_, max_data_value_);
        //min_data_value_ = -10;
        //max_data_value_ = 20;
        float* temp_data = (float*)data_t;
        float value_range = max_data_value_ - min_data_value_;
        cout << "Value Range: " << min_data_value_ << " -> " << min_data_value_ + value_range << endl;
        for ( int i = 0; i < data_size; ++i ) temp_data[i] = ((temp_data[i] - min_data_value_) / value_range);
        break;
                  }
    default:
        OutputErrorLogAndStopRendering("Unsupported data type!");
        break;
    }

    sizes_[0] = sizes_t[0];
    sizes_[1] = sizes_t[1];
    sizes_[2] = sizes_t[2];

    spacings_[0] = spacings_t[0];
    spacings_[1] = spacings_t[1];
    spacings_[2] = spacings_t[2];

    data_format_ = data_format_t;

    //if ( data_ != NULL ) delete data_;
    data_ = data_t;

    UpdateVolumeTexture();

    this->SetViewWindow((min_data_value_ + max_data_value_) / 2, max_data_value_ - min_data_value_);
}

void VolumeRenderWidget::SetBackgroundData(int w, int h, std::vector< float > image_data){
    background_image_data_.assign(image_data.begin(), image_data.end());

    glBindTexture(GL_TEXTURE_2D, background_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_FLOAT, background_image_data_.data());
}

void VolumeRenderWidget::SetViewWindow(float win_center_t, float win_width_t){
    float max_value = win_center_t + win_width_t / 2;
    float min_value = win_center_t - win_width_t / 2;
    if ( max_value > max_data_value_ ) max_value = max_data_value_;
    if ( min_value < min_data_value_ ) min_value = min_data_value_;

    switch (data_format_){
    case GL_UNSIGNED_BYTE:
        volume_bias_ = -1 * (min_value - 0)  / (UCHAR_MAX - 0);
        volume_scale_ = (UCHAR_MAX - 0) / (max_value - min_value);
        break;
    case GL_SHORT:
        volume_bias_ = -1 * (min_value - SHRT_MIN) / (SHRT_MAX - SHRT_MIN) * 2 + 1;
        volume_scale_ = (SHRT_MAX - SHRT_MIN) / (2 * (max_value - min_value));
        break;
    case GL_FLOAT:
        volume_bias_ = -1 * (min_value - min_data_value_) / (max_data_value_ - min_data_value_);
        volume_scale_ = (max_data_value_ - min_data_value_) / (max_value - min_value);
        break;
    default: 
        break;
    }

    this->updateGL();
}

void VolumeRenderWidget::SetSampleStep(float sample_step_t){
    sample_step_ = sample_step_t;

    this->updateGL();
}

QImage* VolumeRenderWidget::GetRenderingImage() {
    this->OptimizeResult();

    QImage* current_image = new QImage(this->width(), this->height(), QImage::Format_ARGB32);
    for (int i = 0; i < this->height(); i++)
        for (int j = 0; j < this->width(); j++) {
            int temp_index = this->width() * (this->height() - 1 - i) + j;
            temp_index *= 4;
            QColor temp_color;
            temp_color.setRedF(optimization_pixel_values_[temp_index]);
            temp_color.setGreenF(optimization_pixel_values_[temp_index + 1]);
            temp_color.setBlueF(optimization_pixel_values_[temp_index + 2]);
            temp_color.setAlphaF(optimization_pixel_values_[temp_index + 3]);
            current_image->setPixel(j, i, temp_color.rgba());
        }

    return current_image;
}

void VolumeRenderWidget::SetTransferFunction(float* values, int entryNumber){
    tf_values_.resize(entryNumber * 4);
    memcpy(tf_values_.data(), values, 4 * entryNumber * sizeof(float));

    UpdateTransferTexture();

    this->updateGL();
}

template < typename T >
void VolumeRenderWidget::GetMinMaxValue(unsigned int data_size, T* data, float& min_value, float& max_value){
    min_value = 1e10;
    max_value = -1e10;
    for ( unsigned int i = 0; i < data_size; ++i ){
        if ( data[i] > max_value ) max_value = data[i];
        if ( data[i] < min_value ) min_value = data[i];
    }
}

void VolumeRenderWidget::OptimizeResult(){
    makeCurrent();
    optimize_run_ = 0;
    is_optimizing_ = true;
    for ( int i = 0; i < 10; ++i ){
        for ( int k = 0; k < random_texture_values_.size(); ++k)
            random_texture_values_[k] = (float)rand() / RAND_MAX;
        UpdateRandomTexture();
        this->paintGL();
    } 
    //is_optimizing_ = false;
}

void VolumeRenderWidget::initializeGL(){
    if ( glewInit() != GLEW_OK ){
        OutputErrorLogAndStopRendering("Initialize OpenGL Failed!");
    }

    glClearColor(1.0, 1.0, 1.0, 1.0);

    InitShaderProgram();
    InitArrayBuffer();
    InitTextures();
    InitFrameBuffer();
}

void VolumeRenderWidget::resizeGL(int w, int h){
    rotation_ball_->SetWindowSize(w, h);

    UpdateBackFaceTexture();
    UpdateOptimizationResultTexture();

    random_texture_values_.resize(this->width() * this->height());
    random_texture_values_.assign(random_texture_values_.size(), 0);
    UpdateRandomTexture();

    optimization_pixel_values_.resize(this->width() * this->height() * 4);
    optimization_pixel_values_.assign(optimization_pixel_values_.size(), 0);

    rendering_result_.resize(this->width() * this->height() * 4);
    rendering_result_.assign(rendering_result_.size(), 0);

    glViewport(0, 0, w, h);
}

void VolumeRenderWidget::paintGL(){
    glClear(GL_COLOR_BUFFER_BIT);
    if ( data_ == NULL ) return;

    if ( background_image_data_.size() != 0 ){
        glUseProgram(render_background_prog_->getId());
        UpdateBackgroundPara();

        glBegin(GL_QUADS);
            glVertex3f(0, 0, 0);
            glVertex3f(1, 0, 0);
            glVertex3f(1, 1, 0);
            glVertex3f(0, 1, 0);
        glEnd();
        if (glGetError() != GL_NO_ERROR ) {
            OutputErrorLogAndStopRendering("Error rendering background!");
        }

        glUseProgram(0);
        glBindVertexArray(0);
    }
    glEnable(GL_CULL_FACE);

    // render the back face texture
    glBindFramebuffer(GL_FRAMEBUFFER, FBO_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, back_face_tex_, 0);

    glUseProgram(back_face_prog_->getId());

    glBindVertexArray(VAO_);
    UpdateBackFacePara();

    glCullFace(GL_FRONT);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, 0);

    if (glGetError() != GL_NO_ERROR ) {
        OutputErrorLogAndStopRendering("Error rendering back face!");
    }
    glBindVertexArray(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glUseProgram(0);

    if ( !is_optimizing_ ){
        glUseProgram(render_volume_prog_->getId());
        glBindVertexArray(VAO_);
        UpdateRenderVolumePara();

        glCullFace(GL_BACK);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, 0);
        if (glGetError() != GL_NO_ERROR ) {
            OutputErrorLogAndStopRendering("Error rendering volume!");
        }

        glUseProgram(0);
        glBindVertexArray(0);
    } else {
        glBindFramebuffer(GL_FRAMEBUFFER, FBO_);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, optimization_result_tex_[0], 0);

        glUseProgram(render_volume_prog_->getId());
        glBindVertexArray(VAO_);
        UpdateRenderVolumePara();

        glCullFace(GL_BACK);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glClear(GL_COLOR_BUFFER_BIT);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, 0);
        if (glGetError() != GL_NO_ERROR ) {
            OutputErrorLogAndStopRendering("Error rendering volume!");
        }

        glUseProgram(0);
        glBindVertexArray(0);

        glReadBuffer(GL_COLOR_ATTACHMENT0);
        glReadPixels(0, 0, this->width(), this->height(), GL_RGBA, GL_FLOAT, rendering_result_.data());

        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        float scale = optimize_run_ / (optimize_run_ + 1);
        for ( int i = 0; i < optimization_pixel_values_.size(); ++i ){
            optimization_pixel_values_[i] += (rendering_result_[i] - optimization_pixel_values_[i]) * (1.0 - scale);
        }
        glDisable(GL_CULL_FACE);
        glViewport(0, 0, this->width(), this->height());
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, 1, 0, 1, 1, 10);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0, 0, -3);

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
        glEnable(GL_BLEND); 

        glActiveTexture(GL_TEXTURE0);
        glEnable(GL_TEXTURE_2D);
        glColor4f(1.0, 1.0, 1.0, 1.0);
        glBindTexture(GL_TEXTURE_2D, optimization_result_tex_[1]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->width(), this->height(), 0, GL_RGBA, GL_FLOAT, optimization_pixel_values_.data());
        glBegin(GL_QUADS);
        glTexCoord2f(0, 0); glVertex3f(0, 0, 0);
        glTexCoord2f(0, 1); glVertex3f(0, 1, 0);
        glTexCoord2f(1, 1); glVertex3f(1, 1, 0);
        glTexCoord2f(1, 0); glVertex3f(1, 0, 0);
        glEnd();
        glDisable(GL_TEXTURE_2D);
        glDisable(GL_BLEND);

        optimize_run_++;
    }

    glDisable(GL_CULL_FACE);
}

void VolumeRenderWidget::UpdateBackgroundPara(){
    const glm::mat4 model = glm::mat4(rotation_ball_->Rotation() * glm::mat3(glm::scale(glm::vec3(zoom_factor_))));
    const glm::mat4 modelView = glm::lookAt(glm::vec3(0.0f, 0.0f, 4.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * model;
    const glm::mat4 modelViewProj = glm::perspective(45.0f, (float)this->width() / this->height(), 0.1f, 10.0f) * modelView;

    GLuint location = glGetUniformLocation(render_background_prog_->getId(), "ModelViewProj");
    if (location != -1) 
        glUniformMatrix4fv(location, 1, GL_FALSE, &modelViewProj[0][0]);

    if (location == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Unable to get location for modelViewProj in back face rendering!");
    }

    float local_size[3];
    float max_size = 0;
    for ( int i = 0; i < 3; ++i ) {
        local_size[i] = sizes_[i] * spacings_[i];
        if ( max_size < local_size[i] ) max_size = local_size[i];
    }
    if ( max_size != 0 ){
        for ( int i = 0; i < 3; ++i ) local_size[i] /= max_size;
    }
    location = glGetUniformLocation(render_background_prog_->getId(), "VertexScale");
    if (location != -1)
        glUniform3fv(location, 1, static_cast<GLfloat *>(local_size));

    if (location == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Unable to get location for VertexScale in back face rendering!");
    }

    location = glGetUniformLocation(render_background_prog_->getId(), "BackgroundTex");
    if ( location != -1 ){
        glUniform1i(location , 0);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, background_tex_);
    }
    GLenum err = glGetError();
    if ( location == -1 || glGetError() != GL_NO_ERROR ){
        OutputErrorLogAndStopRendering("Error binding texture in rendering background!");
    }
}

void VolumeRenderWidget::UpdateBackFacePara(){
    const glm::mat4 model = glm::mat4(rotation_ball_->Rotation() * glm::mat3(glm::scale(glm::vec3(zoom_factor_))));
    const glm::mat4 modelView = glm::lookAt(glm::vec3(0.0f, 0.0f, 4.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * model;
    const glm::mat4 modelViewProj = glm::perspective(45.0f, (float)this->width() / this->height(), 0.1f, 10.0f) * modelView;
    
    GLuint location = glGetUniformLocation(back_face_prog_->getId(), "ModelViewProj");
    if (location != -1) 
        glUniformMatrix4fv(location, 1, GL_FALSE, &modelViewProj[0][0]);

    if (location == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Unable to get location for modelViewProj in back face rendering!");
    }

    float local_size[3];
    float max_size = 0;
    for ( int i = 0; i < 3; ++i ) {
        local_size[i] = sizes_[i] * spacings_[i];
        if ( max_size < local_size[i] ) max_size = local_size[i];
    }
    if ( max_size != 0 ){
        for ( int i = 0; i < 3; ++i ) local_size[i] /= max_size;
    }
    location = glGetUniformLocation(back_face_prog_->getId(), "VertexScale");
    if (location != -1)
        glUniform3fv(location, 1, static_cast<GLfloat *>(local_size));

    if (location == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Unable to get location for VertexScale in back face rendering!");
    }

    location = glGetAttribLocation(back_face_prog_->getId(), "PositionIn");
    if (location != -1){
        glBindBuffer(GL_ARRAY_BUFFER, VBO_);
        glVertexAttribPointer(location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(location);
    }
}

void VolumeRenderWidget::UpdateRenderVolumePara(){
    const glm::mat4 model = glm::mat4(rotation_ball_->Rotation() * glm::mat3(glm::scale(glm::vec3(zoom_factor_))));
    const glm::mat4 modelView = glm::lookAt(glm::vec3(0.0f, 0.0f, 4.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * model;
    const glm::mat4 modelViewProj = glm::perspective(45.0f, (float)this->width() / this->height(), 0.1f, 10.0f) * modelView;
    glm::mat3 upper_left_modelView;
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j ) 
            upper_left_modelView[i][j] = modelView[i][j];
    const glm::mat3 normalMatrix = glm::transpose(glm::inverse(upper_left_modelView));

    light_half_vec_ = glm::normalize((glm::vec3(0, 0, 1) - light_dir_)); 

    GLuint location = glGetUniformLocation(render_volume_prog_->getId(), "ModelViewProj");
    if (location != -1) 
        glUniformMatrix4fv(location, 1, GL_FALSE, &modelViewProj[0][0]);
    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating ModelViewProj in rendering volume!");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "NormalMatrix");
    if (location != -1) 
        glUniformMatrix3fv(location, 1, GL_FALSE, &normalMatrix[0][0]);

    if (location == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Unable to get location for normal matrix in volume rendering!");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "TransferFunc");
    if ( location != -1 ){
        glUniform1i(location , 0);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_1D, transfer_function_tex_);
    }
    if ( location == -1 || glGetError() != GL_NO_ERROR ){
        OutputErrorLogAndStopRendering("Error binding texture in rendering volume!");
    }
    location = glGetUniformLocation(render_volume_prog_->getId(), "BackFaceTex");
    if ( location != -1 ){
        glProgramUniform1i(render_volume_prog_->getId(), location , 1);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, back_face_tex_);
    }
    if ( location == -1 || glGetError() != GL_NO_ERROR ){
        OutputErrorLogAndStopRendering("Error binding texture in rendering volume!");
    }
    location = glGetUniformLocation(render_volume_prog_->getId(), "VolumeTex");
    if ( location != -1 ){
        glProgramUniform1i(render_volume_prog_->getId(), location , 2);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_3D, volume_tex_);
    }
    if ( location == -1 || glGetError() != GL_NO_ERROR ){
        OutputErrorLogAndStopRendering("Error binding texture in rendering volume!");
    }
    location = glGetUniformLocation(render_volume_prog_->getId(), "RandomTex");
    if ( location != -1 ){
        glProgramUniform1i(render_volume_prog_->getId(), location , 3);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, random_tex_);
    }
    if ( location == -1 || glGetError() != GL_NO_ERROR ){
        OutputErrorLogAndStopRendering("Error binding texture in rendering volume!");
    }

    float local_size[3];
    float max_size = 0;
    for ( int i = 0; i < 3; ++i ) {
        local_size[i] = sizes_[i] * spacings_[i];
        if ( max_size < local_size[i] ) max_size = local_size[i];
    }
    if ( max_size != 0 ){
        for ( int i = 0; i < 3; ++i ) local_size[i] /= max_size;
    }
    location = glGetUniformLocation(render_volume_prog_->getId(), "VertexScale");
    if (location != -1)
        glUniform3fv(location, 1, static_cast<GLfloat *>(local_size));
    if (location == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error updating voxel size!");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "StepSize");
    if (location != -1)
        glUniform1f(location, sample_step_);
    if (location == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error updating sampling step size");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "VolumeScale");
    if (location != -1)
        glUniform1f(location, volume_scale_);
    if (location == -1 || glGetError() != GL_NO_ERROR)
    {
        OutputErrorLogAndStopRendering("Error updating volume scale!");
        return;
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "VolumeBias");
    if (location != -1)
        glUniform1f(location, volume_bias_);
    if (location == -1 || glGetError() != GL_NO_ERROR)
    {
        OutputErrorLogAndStopRendering("Error updating volume bias");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "ScreenResolution");
    if (location != -1){
        float resolution[2];
        resolution[0] = this->width();
        resolution[1] = this->height();
        glUniform2fv(location, 1, resolution);
    }

    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating screen resolution!");
    }

    /*location = glGetUniformLocation(render_volume_prog_->getId(), "BackgroundColor");
    if (location != -1){
        glUniform3fv(location, 1, &back_ground_color_[0]);
    }

    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating background color!");
    }*/

    location = glGetUniformLocation(render_volume_prog_->getId(), "LightDir");
    if (location != -1){
        glUniform3fv(location, 1, &light_dir_[0]);
    }

    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating LightDir!");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "LightDiffuse");
    if (location != -1){
        glUniform4fv(location, 1, &light_diffuse_[0]);
    }

    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating LightDiffuse!");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "LightAmbient");
    if (location != -1){
        glUniform4fv(location, 1, &light_ambient_[0]);
    }

    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating LightAmbient!");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "LightSpecular");
    if (location != -1){
        glUniform4fv(location, 1, &light_specular_[0]);
    }

    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating LightSpecular!");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "Shininess");
    if (location != -1){
        glUniform1f(location, shininess_);
    }

    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating Shininess!");
    }

    location = glGetUniformLocation(render_volume_prog_->getId(), "LightHalfVec");
    if (location != -1){
        glUniform3fv(location, 1, &light_half_vec_[0]);
    }

    if (location == -1 || glGetError() != GL_NO_ERROR) {
        OutputErrorLogAndStopRendering("Error updating LightSpecular!");
    }

    location = glGetAttribLocation(render_volume_prog_->getId(), "PositionIn");
    if (location != -1){
        glBindBuffer(GL_ARRAY_BUFFER, VBO_);
        glVertexAttribPointer(location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(location);
    }
}

void VolumeRenderWidget::InitShaderProgram(){
    GlShader* back_vert_shader = new GlShader("./Shaders/back_position.vsh");
    GlShader* back_frag_shader = new GlShader("./Shaders/back_position.psh");
    if (back_vert_shader->status() != GlShader::OK || back_frag_shader->status() != GlShader::OK || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Building back face shader error!");
    }
    back_face_prog_ = new GlProgram(back_vert_shader, back_frag_shader);
    if (back_face_prog_->buildProgram() == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering(back_face_prog_->getInfoLog()->data());
    }

    GlShader* volume_vert_shader = new GlShader("./Shaders/volume_rendering.vsh");
    GlShader* volume_frag_shader = new GlShader("./Shaders/volume_rendering_with_light.psh");
    if (volume_vert_shader->status() != GlShader::OK || volume_frag_shader->status() != GlShader::OK || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Building volume shader error!");
    }
    render_volume_prog_ = new GlProgram(volume_vert_shader, volume_frag_shader);
    if (render_volume_prog_->buildProgram() == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering(render_volume_prog_->getInfoLog()->data());
    }

    GlShader* background_vert_shader = new GlShader("./Shaders/background.vsh");
    GlShader* background_frag_shader = new GlShader("./Shaders/background.psh");
    if (background_vert_shader->status() != GlShader::OK || background_frag_shader->status() != GlShader::OK || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Building volume shader error!");
    }
    render_background_prog_ = new GlProgram(background_vert_shader, background_frag_shader);
    if (render_background_prog_->buildProgram() == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering(render_background_prog_->getInfoLog()->data());
    }
}

void VolumeRenderWidget::InitArrayBuffer(){
    const glm::vec3 vertices[8] = {
        glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 1.0f),
        glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 1.0f, 1.0f), glm::vec3(1.0f, 1.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f)
    };

    glGenBuffers(1, &VBO_);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    if (VBO_ == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error during generating array buffer!");
    }

    const GLushort indices[36] = {
        2, 3, 0, 0, 3, 1, // bottom
        0, 1, 4, 4, 1, 5, // left
        6, 2, 4, 4, 2, 0, // rear
        6, 7, 2, 2, 7, 3, // right
        3, 7, 1, 1, 7, 5, // front
        4, 5, 6, 6, 5, 7, // top
    };

    glGenBuffers(1, &EBO_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    if (EBO_ == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Element Array Buffer Generation Error");
    }

    glGenVertexArrays(1, &VAO_);
    glBindVertexArray(VAO_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_);

    if (VAO_ == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Vertex Array Buffer Generation Error");
    }

    glBindVertexArray(0);
}

void VolumeRenderWidget::InitTextures(){
    glGenTextures(1, &transfer_function_tex_);
    glBindTexture(GL_TEXTURE_1D, transfer_function_tex_);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 0, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_1D, 0);

    if (glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error when initializing texture for transfer function!");
    }

    glGenTextures(1, &back_face_tex_);
    glBindTexture(GL_TEXTURE_2D, back_face_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error when initializing texture for back space!");
    }

    glGenTextures(1, &background_tex_);
    glBindTexture(GL_TEXTURE_2D, background_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error when initializing texture for background!");
    }

    glGenTextures(1, &random_tex_);
    glBindTexture(GL_TEXTURE_2D, random_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, 0, 0, 0, GL_ALPHA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error when random texture!");
    }

    glGenTextures(2, optimization_result_tex_);
    glBindTexture(GL_TEXTURE_2D, optimization_result_tex_[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glBindTexture(GL_TEXTURE_2D, optimization_result_tex_[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error when initializing texture for back space!");
    }

    glGenTextures(1, &volume_tex_);
    glBindTexture(GL_TEXTURE_3D, volume_tex_);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_ALPHA, 0, 0, 0, 0, GL_ALPHA, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_3D, 0);

    if (glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error when initializing texture for volume!");
    }
}

void VolumeRenderWidget::InitFrameBuffer(){
    glGenFramebuffers(1, &FBO_);

    if (FBO_ == -1 || glGetError() != GL_NO_ERROR){
        OutputErrorLogAndStopRendering("Error when generating frame buffer!");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, FBO_);

    if ( glGetError() != GL_NO_ERROR ){
        OutputErrorLogAndStopRendering("Error when binding color attachment!");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void VolumeRenderWidget::UpdateVolumeTexture(){
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glBindTexture(GL_TEXTURE_3D, volume_tex_);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_ALPHA, sizes_[0], sizes_[1], sizes_[2], 0, GL_ALPHA, data_format_, data_);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_3D, 0);
}

void VolumeRenderWidget::UpdateTransferTexture(){
    glBindTexture(GL_TEXTURE_1D, transfer_function_tex_);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, tf_values_.size() / 4, 0, GL_RGBA, GL_FLOAT, tf_values_.data());
    glBindTexture(GL_TEXTURE_1D, 0);
}

void VolumeRenderWidget::UpdateBackFaceTexture(){
    glBindTexture(GL_TEXTURE_2D, back_face_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->width(), this->height(), 0, GL_RGBA, GL_FLOAT, NULL);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void VolumeRenderWidget::UpdateOptimizationResultTexture(){
    glBindTexture(GL_TEXTURE_2D, optimization_result_tex_[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->width(), this->height(), 0, GL_RGBA, GL_FLOAT, NULL);
    glBindTexture(GL_TEXTURE_2D, optimization_result_tex_[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->width(), this->height(), 0, GL_RGBA, GL_FLOAT, NULL);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void VolumeRenderWidget::UpdateRandomTexture(){
    glBindTexture(GL_TEXTURE_2D, random_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, this->width(), this->height(), 0, GL_ALPHA, GL_FLOAT, random_texture_values_.data());
    glBindTexture(GL_TEXTURE_2D, 0);
}

void VolumeRenderWidget::OutputErrorLogAndStopRendering(const char* message){
    QMessageBox::information(this, tr("Warning"), QString::fromLocal8Bit(message));
    exit(-1);
}

void VolumeRenderWidget::mousePressEvent(QMouseEvent *event){
    rotation_ball_->Begin(event->x(), event->y());
}

void VolumeRenderWidget::mouseMoveEvent(QMouseEvent *event){
    rotation_ball_->Drag(event->x(), event->y());

    RestartOptimizingTimer();
}

void VolumeRenderWidget::mouseReleaseEvent(QMouseEvent *event){

}

void VolumeRenderWidget::wheelEvent(QWheelEvent *event){
    zoom_factor_ += event->delta() / 1000.0f;

    if ( zoom_factor_ < 0.1 ) zoom_factor_ = 0.1;
    if ( zoom_factor_ > 10 ) zoom_factor_ = 10;

    RestartOptimizingTimer();
}

void VolumeRenderWidget::OnOptimizingTimerTimeOut(){
    for ( int k = 0; k < random_texture_values_.size(); ++k)
        random_texture_values_[k] = (float)rand() / RAND_MAX;
    UpdateRandomTexture();
    this->updateGL();
    if ( optimize_run_ < 10 ) optimizing_timer_->start(200);
}

void VolumeRenderWidget::RestartOptimizingTimer(){
    is_optimizing_ = true;
    optimize_run_ = 0;
    optimization_pixel_values_.assign(optimization_pixel_values_.size(), 0);

    OnOptimizingTimerTimeOut();
}