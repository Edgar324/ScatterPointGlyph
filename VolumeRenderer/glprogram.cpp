#include "glprogram.h"
#include "glshader.h"

GlProgram::GlProgram(GlShader* vertexShader, GlShader* fragmentShader)
{
	programId_ = -1;
	infoLog_ = new std::vector< GLchar >;
	infoLog_->clear();

	vertexShader_ = vertexShader;
	fragmentShader_ = fragmentShader;

	status_ = INITIALIZED;
}

GlProgram::~GlProgram()
{
	if (vertexShader_ != 0) delete vertexShader_;
	if (fragmentShader_ != 0) delete fragmentShader_;

	infoLog_->clear();
	delete infoLog_;

	if (programId_ != -1) glDeleteProgram(programId_);
}

GLuint GlProgram::buildProgram()
{
	programId_ = glCreateProgram();

	const GLchar* vertexShaderSource = vertexShader_->getContent();
	GLuint vertexId = compileShader(GL_VERTEX_SHADER, &vertexShaderSource);
	if (vertexId != -1)
	{
		glAttachShader(programId_, vertexId);
		glDeleteShader(vertexId);
	}
	else
	{
		status_ = VERTEX_SHADER_ERROR;
		return -1;
	}
	
	const GLchar* fragmentShaderSource = fragmentShader_->getContent();
	GLuint fragmentId = compileShader(GL_FRAGMENT_SHADER, &fragmentShaderSource);
	if (fragmentId != -1)
	{
		glAttachShader(programId_, fragmentId);
		glDeleteShader(fragmentId);
	}
	else
	{
		status_ = FRAGMENT_SHADER_ERROR;
		return -1;
	}

	glLinkProgram(programId_);
	if (!checkProgramStatus(programId_))
	{
		glDeleteShader(programId_);
		programId_ = -1;
		status_ = LINK_ERROR;
		return -1;
	}
	
	status_ = OK;

	return programId_;
}

GLuint GlProgram::compileShader(GLenum shaderType, const GLchar** shaderSource)
{
	const GLuint name = glCreateShader(shaderType);

	glShaderSource(name, 1, shaderSource, 0);
	glCompileShader(name);

	if (!checkShaderStatus(name))
	{
		glDeleteShader(name);
		return -1;
	}
	else
		return name;
	
}

bool GlProgram::checkShaderStatus(GLuint id)
{
	GLint status;
	glGetShaderiv(id, GL_COMPILE_STATUS, &status);

	if (!status) 
	{
		GLint s;
		glGetShaderiv(id, GL_INFO_LOG_LENGTH, &s);
		infoLog_->resize(s);
		glGetShaderInfoLog(id, static_cast<GLsizei>(infoLog_->size()), NULL, &infoLog_->at(0));

		return false;
	}

	return true;
}

bool GlProgram::checkProgramStatus(GLuint id)
{
	GLint status;
	glGetProgramiv(id, GL_LINK_STATUS, &status);

	if (!status) 
	{
		GLint s;
		glGetProgramiv(id, GL_INFO_LOG_LENGTH, &s);
		infoLog_->resize(s);
		glGetProgramInfoLog(id, static_cast<GLsizei>(infoLog_->size()), NULL, &infoLog_->at(0));
		return false;
	}
	else
		return true;
}