#pragma once

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GLM/glm.hpp>
#include <GLM/vec3.hpp>
#include <GLM/vec4.hpp>
#include <GLM/mat4x4.hpp>
#include <GLM/gtc/matrix_transform.hpp>
#include <GLM/gtx/rotate_vector.hpp>
#include <GLM/gtc/type_ptr.hpp>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cfloat>

#define PI 3.14159f
typedef unsigned char uchar8;



//material

struct Material
{
	Material(glm::vec3 amb, glm::vec3 diff, glm::vec3 spec, float shin) : ambient(amb), diffuse(diff), specular(spec), shininess(shin) {};
	Material(const Material& _m) { ambient = _m.ambient; diffuse = _m.diffuse; specular = _m.specular; shininess = _m.shininess; };
	glm::vec3 ambient;
	glm::vec3 diffuse;
	glm::vec3 specular;
	float shininess;
};

//Nome correto deveria ter sido Point em vez de vertex2D
struct Vertex2D
{
	Vertex2D() : x(0.0), y(0.0) {};
	Vertex2D(float _x, float _y) : x(_x), y(_y) {};
	float x;
	float y;
};


//struct vertex customizado para as primeiras iterações do programa
//TODO: implementar UV , smooth de normais e indices
struct Vertex
{
	Vertex(){};
	
	Vertex(const glm::vec3& _position, const glm::vec3& _normal)
		: position(_position)
		, normal(_normal)
		, uv(0)
	{};

	Vertex(const glm::vec3& _position, const glm::vec3& _normal, const glm::vec2& _uv)
		: position(_position)
		, normal(_normal)
		,uv(_uv) {};

	Vertex(const glm::vec3& _position,const glm::vec2& _uv)
		: position(_position)
		, normal(0)
		, uv(_uv) {};

	Vertex(const Vertex& _v)
		:position(_v.position)
		,normal(_v.normal)
		,uv(_v.uv)
	{};


	glm::vec3 position;
	glm::vec3 normal;
	glm::vec2 uv;
};






