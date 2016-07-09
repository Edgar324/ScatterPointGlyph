#ifndef GLPROGRAM_H
#define GLPROGRAM_H

#include "GL/glew.h"
#include <vector>

class GlShader;

class GlProgram
{
public:
	GlProgram(GlShader* vertexShader, GlShader* fragmentShader);
	~GlProgram();

	enum ProgramStatus{INITIALIZED = 0, VERTEX_SHADER_ERROR, GEOMETRY_SHADER_ERROR, FRAGMENT_SHADER_ERROR, LINK_ERROR, OK};

	GLuint buildProgram();
	GLuint getId() { return programId_; }
	ProgramStatus status() { return status_; }
	std::vector< GLchar >* getInfoLog() { return infoLog_; }

private:
	GlShader* vertexShader_;
	GlShader* fragmentShader_;

	GLuint programId_;

	ProgramStatus status_;
	std::vector< GLchar >* infoLog_;

	GLuint compileShader(GLenum shaderType, const GLchar** shaderSource);
	bool checkShaderStatus(GLuint id);
	bool checkProgramStatus(GLuint id);
};

#endif