#ifndef GLSHADER_H
#define GLSHADER_H

#include "GL/glew.h"

class GlShader
{
public:
	GlShader(const char*);
	~GlShader();

	enum ShaderStatus {INITIALIZED = 0, FILE_NOT_LOADED, MEMORY_MALLOC_ERROR, OK};

	ShaderStatus status() { return shaderStatus_; }

	GLchar* getContent() { return content_; }

private:
	ShaderStatus shaderStatus_;
	GLchar* content_;
};

#endif