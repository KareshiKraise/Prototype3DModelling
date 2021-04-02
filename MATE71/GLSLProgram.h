#ifndef _SHADERPROG_
#define _SHADERPROG_

#include <string>
#include "Includes.h"



	class GLSLProgram
	{
	public:
		GLSLProgram();
		~GLSLProgram();

		void compileShaders(const std::string& vertexShaderFilePath, const std::string& fragmentShaderFilePath);
		void compileShadersFromSource(const char* vertexSource, const char* fragmentSource);

		void linkShaders();
		void addAttribute(const char* attributeName);

		void use();
		void unuse();

		GLint getUniformLocation(const char* uniformName);

		void dispose();


	private:

		int _numAttributes;
		
		void compileShader(const char* source, const std::string& name, GLuint id);
		


		GLuint _programID;
		GLuint _vertexShaderID;
		GLuint _fragmentShaderID;
		







	};

#endif