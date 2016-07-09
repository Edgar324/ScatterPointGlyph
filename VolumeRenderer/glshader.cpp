#include "glshader.h"
#include <fstream>

GlShader::GlShader(const char* fileName)
{
	shaderStatus_ = INITIALIZED;

	std::ifstream file(fileName, std::ios::binary);
	if (!file)
	{
		shaderStatus_ = FILE_NOT_LOADED;
	}

	file.seekg(0, std::ios::end);
	std::streamoff length = file.tellg();
	file.seekg(0, std::ios::beg);

	content_ = new GLchar[static_cast< unsigned int >(length + 1)];
	if (content_ == NULL)
	{
		shaderStatus_ = MEMORY_MALLOC_ERROR;
	}
	else
	{
		file.read(content_, length);
		content_[length] = 0;
	}

	file.close();

	if (shaderStatus_ == INITIALIZED) shaderStatus_ = OK;
}

GlShader::~GlShader()
{
	if (content_ != NULL) delete content_;
}