#define _CRT_SECURE_NO_WARNINGS

//TODO: codigo sujo e necessitando de reorganizacao urgente


#include <iostream>
#include "Includes.h"
#include "GLSLProgram.h"
#include <SOIL/SOIL.h>
/*-------------------------------------------------------------------------------------------------------------------------*/
//GLOBAIS
float WIDTH  = 1200;
float HEIGHT = 600;
bool is_parametric = true;

int has_granite = 0;
int has_metal = 0;

static glm::mat4 Model = glm::mat4(1.0);

Material metal(glm::vec3(0.04, 0.04, 0.04), glm::vec3(0.75, 0.75, 0.75), glm::vec3(0.9, 0.9, 0.9), 32.0f);

Material current(metal);

float FOV = 45;
int last_mx = 0, last_my = 0, cur_mx = 0, cur_my = 0;
int arcball_on = false;
glm::vec3 eyePos = glm::vec3(150.0f, 0.0f, 0.0f);
glm::vec3 lightPos = glm::vec3(150.0f, 0.0f, 0.0f);
glm::vec3 objectColor = glm::vec3(0.5, 0.5, 0.5);
glm::vec3 lightColor = glm::vec3(1.0, 1.0, 1.0);

float specular = 0.95f;

bool see_X = false;
bool see_Y = false;
bool see_Z = false;
bool see_Iso = false;


bool carry_point = false;
unsigned int c_index;

struct Boundaries
{	
	glm::vec3 min;
	glm::vec3 max;
	glm::vec3 center;

} Bounds;

static glm::mat4 Ortho = glm::ortho(-300.0f, 300.f, -300.0f, 300.0f, 300.0f, -300.0f);

static glm::mat4 View = glm::lookAt(eyePos,
									glm::vec3(0.0f, 0.0f, 0.0f),
									glm::vec3(0.0f, 1.0f, 0.0f));

static glm::mat4 Projection = glm::perspective(FOV, (WIDTH / 2) / (HEIGHT), 300.0f, -300.0f);


//glm::vec2 comp;
glm::mat3 last_rotation = glm::mat3(1.0);
glm::mat3 cur_rotation = glm::mat3(1.0);

GLuint axisVBO;
GLuint axisVAO;

GLuint pointsVAO;
GLuint pointsVBO;

GLuint curveVAO;
GLuint curveVBO;

GLuint solidVAO;
GLuint solidVBO;

GLuint axis3DVAO;
GLuint axis3DVBO;

GLuint cpointVAO;
GLuint cpointVBO;

GLint matrixID;
GLint uniID;

GLSLProgram AXIS_PROGRAM;
GLSLProgram POINTS_PROGRAM;
GLSLProgram CURVE_PROGRAM;
GLSLProgram SOLID_PROGRAM;
GLSLProgram AXIS3D_PROGRAM;


std::vector<glm::vec2> control_points; //control points
std::vector<glm::vec2> p_curve; //parametric bezier curve
std::vector<glm::vec3> lathe_coords;

std::vector<Vertex> model;

GLuint texID;

//int texSize = 128;
//GLuint texName = 0;
//GLubyte *texPtr;

bool fullscreen3D = false;
bool fullscreen2D = false;

//eixo 3d 
static const GLfloat axis_data2[] =
{
	-WIDTH / 4,    0.0f,    0.0f, //begin horizontal line
	WIDTH / 4,    0.0f,    0.0f,  //end horizontal line

	0.0f , -HEIGHT / 2,    0.0f,  //begin vertical line
	0.0f ,  HEIGHT / 2,    0.0f,  //end vertical line

	0.0f ,    0.0f, -300.0f,	  //z axis line begin	
	0.0f ,    0.0f,  300.0f,	  //z axis line end

	1.0f,    0.0f,    0.0f,       //color red
	1.0f,    0.0f,   0.0f,		  //color red	

	0.0f ,    1.0f,    0.0f,      //color green
	0.0f ,    1.0f,    0.0f,      //color green

	0.0f ,    0.0f,    1.0f,      //color blue
	0.0f ,    0.0f,    1.0f     //color blue
	
};

//eixo 2d
static GLfloat axis_data[] =
{
	-WIDTH / 4 ,    0.0f,    0.0f,       //begin horizontal line
	WIDTH / 4 ,    0.0f,    0.0f,	    //end horizontal line
	0.0f , -HEIGHT / 2,    0.0f,       //begin vertical line
	0.0f ,  HEIGHT / 2,    0.0f,       //end vertical line
	0.0f ,    0.0f, -300.0f,		//z axis line begin	
	0.0f ,    0.0f,  300.0f		//z axis line end
};



//Para solido de revolucao
const int   SECTIONS = 45;
const double STEP = (PI * 2) / (double)SECTIONS;

//Para solido parametrico
const int SEGMENT = 12;
const double ANGLE = (PI * 2) / (double)SEGMENT;


std::vector<std::vector<glm::vec3>> rotatedCpoints;

//FIM DAS GLOBAIS
/*-------------------------------------------------------------------------------------------------------------------------*/


//Calcula aproximacao das UV's usando coordenadas esfericas
void sphericalCoords(std::vector<Vertex>& _vert)
{
	double fi;
	double alpha;
	
	
	//std::vector<glm::vec2> list_spherical;
	std::vector<glm::vec2> my_uv;

	glm::vec3 unit;
	for (int i = 0; i < _vert.size(); i++)
	{
		unit = glm::normalize(_vert[i].position);

		//Eu sinceramente nao entendo porque isso está funcionando com (PI), as fontes dizem que o angulo alpha deveria ser divido por 2*PI e que atan2 eh uma opcao melhor
		//Fiz tweaks adhoc mas nao estou satisfeito
				
		fi = 0.5 - (asin(unit.y) / PI);
		alpha = 0.5 + (atan(unit.z / unit.x)) / PI;
		

		_vert[i].uv = (glm::vec2(alpha, fi));
	}
	
		
}

//Gera um solido de revolução
std::vector<Vertex> genSOR(const std::vector<glm::vec2>& points2D)
{
	std::vector<glm::vec3> aux(points2D.size());

	for (int i = 0; i < points2D.size(); i++)
	{
		aux[i] = glm::vec3(points2D[i].x, points2D[i].y, 0);
	}

	std::vector<std::vector<glm::vec3>> layers(SECTIONS, std::vector<glm::vec3>(aux.size()));

	for (int i = 0; i < layers.size(); i++)
	{
		for (int j = 0; j < aux.size(); j++)
		{
			layers[i][j] = glm::vec3(aux[j].x * cos(STEP * i) + aux[j].z * sin(STEP * i),
									 aux[j].y,
									 aux[j].z * cos(STEP * i) - aux[j].x * sin(STEP * i));
		}	
	}
	layers.shrink_to_fit();

	std::vector<Vertex> verts;

	//std::vector<Vertex> test;

	glm::vec2 uv = glm::vec2(0);
	float uv_step = 1 / SECTIONS;

	for (int i = 1; i < layers.size(); ++i)
	{
		const std::vector<glm::vec3>& prvLayer = layers[i - 1];
		const std::vector<glm::vec3>& curLayer = layers[i - 0];

		for (int j = 1; j < aux.size(); ++j)
		{
			const glm::vec3& LL = prvLayer[j - 1];
			const glm::vec3& LR = prvLayer[j - 0];
			const glm::vec3& UL = curLayer[j - 1];
			const glm::vec3& UR = curLayer[j - 0];


			const glm::vec3 normal0 = glm::normalize(glm::cross(UR - LL, UL - LL));
			
			verts.push_back(Vertex(LL, normal0));
			verts.push_back(Vertex(UR, normal0));
			verts.push_back(Vertex(UL, normal0));


			const glm::vec3 normal1 = glm::normalize(glm::cross(LR - LL, UL - LL));
			verts.push_back(Vertex(LL, normal1));
			verts.push_back(Vertex(LR, normal1));
			verts.push_back(Vertex(UR, normal1));
						

		}
	}


	const std::vector<glm::vec3>& prvLayer = layers[layers.size() - 1];
	const std::vector<glm::vec3>& curLayer = layers[0];
	for (int j = 1; j < aux.size(); ++j)
	{
		const glm::vec3& LL = prvLayer[j - 1];
		const glm::vec3& LR = prvLayer[j - 0];
		const glm::vec3& UL = curLayer[j - 1];
		const glm::vec3& UR = curLayer[j - 0];


		const glm::vec3 normal0 = glm::normalize(glm::cross(UR - LL, UL - LL));
		verts.push_back(Vertex(LL, normal0));
		verts.push_back(Vertex(UR, normal0));
		verts.push_back(Vertex(UL, normal0));


		const glm::vec3 normal1 = glm::normalize(glm::cross(LR - LL, UL - LL));
		verts.push_back(Vertex(LL, normal1));
		verts.push_back(Vertex(LR, normal1));
		verts.push_back(Vertex(UR, normal1));

	}


	// flip normals
	for (int i = 0; i < verts.size(); i++)
	{
		verts[i].normal = (glm::dot(verts[i].normal, verts[i].position) < 0.0) ? verts[i].normal : -verts[i].normal;
	}


	return verts;

}

//Calculo de !(fatorial)
unsigned int fact(unsigned int n)
{
	if (n <= 1)
	{
		return 1;
	}

	unsigned int fac = 1;
	for (int i = 1; i <= n; i++)
	{
		fac *= i;
	}
	return fac;
}

//calcula component bernstein para curva
double bern_poly(double a, double b)
{
	return fact(a) / (fact(b) * fact(a - b));
}



/*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

//Binomio de newton para superficie
float binomial_coefficient(int i, int n)
{
	if (i >= 0 && n >=0)
	{
		return 1.0f * fact(n) / (fact(i) * fact(n - i));
	}
	
	return 1;
}

//calcula component bernstein para superficie
double bern_poly2(int i, int n, float u)
{
	return binomial_coefficient(i, n) * powf(u, i) * powf(1 - u, n - i);
}

//Calcula uma curva de bezier
glm::vec2 evaluate_curve(const std::vector<glm::vec2>& c_points, float t, int degree)
{
	glm::vec2 bez_P;	//bezier 
	int n = degree;
	
	for (int i = 0; i <= n; i++)
	{
		bez_P.x += c_points[i].x * bern_poly(n, i) * pow(1-t,n-i) * pow(t,i);
		bez_P.y += c_points[i].y * bern_poly(n, i) * pow(1 - t, n - i) * pow(t, i);
	}
	
	return bez_P;
}

//Rotaciona somente os control points da superficie parametrica
std::vector<std::vector<glm::vec3>> rotate_points(const std::vector<glm::vec2>& points2D)
{
	std::vector<glm::vec3> aux(points2D.size());

	for (int i = 0; i < points2D.size(); i++)
	{
		aux[i] = glm::vec3(points2D[i].x, points2D[i].y, 0);
	}

	std::vector<std::vector<glm::vec3>> slices(SEGMENT, std::vector<glm::vec3>(aux.size()));

	
	slices[0] = aux;
		
	for (int j = 1; j < (slices.size()); j++)
	{
		for (int k = 0; k < slices[j].size(); k++)
		{
			//slices[j][k] = glm::vec3(aux[k].z * sin(ANGLE * j) + aux[k].x * cos(ANGLE * j),
			//						 aux[k].y,
			//				         aux[k].z * cos(ANGLE * j) - aux[k].x * sin(ANGLE * j));

			slices[j][k] = glm::vec3(slices[j-1][k].z * sin(ANGLE) + slices[j-1][k].x * cos(ANGLE),
								     slices[j-1][k].y,
									 slices[j-1][k].z * cos(ANGLE) - slices[j-1][k].x * sin(ANGLE));

		}
	}


	



	
	slices[SEGMENT - 1] = slices[0];

	//slices.shrink_to_fit();
	//for (int i = 0; i < slices.size(); i++)
	//{
	//	slices[i].shrink_to_fit();
	//}
	
	/*for (int j = 0; j < slices.size(); j++)
	{
		for (int k = 0; k < slices[j].size(); k++)
		{
			std::cout << slices[j][k].x << " " << slices[j][k].y << " " << slices[j][k].z << std::endl;
		}
	}*/
	
	return slices;

}

//calcula uma superficie de bezier
std::vector <Vertex> bez_surface(const std::vector<std::vector<glm::vec3>>& c_points)
{
	int M = c_points.size() - 1;
	int N = c_points[0].size() - 1;

	//Resolucoes maiores que M*10 / N*10 causam lag
	int M_res = M * 4;
	int N_res = N * 4;

	
	double u = 0;
	double v = 0;

	std::vector<std::vector<Vertex>> out(M_res, std::vector<Vertex>(N_res));

	for (int i = 0; i < M_res; i++)
	{
		u = (double)i / (double)(M_res - 1);

		for (int j = 0; j < N_res; j++)
		{
			v = (double)j / (double)(N_res - 1);
			
			for (int ki = 0; ki <= M; ki++)
			{
				double bi = bern_poly2(ki, M, u);

				for (int kj = 0; kj <= N; kj++)
				{
					double bj = bern_poly2(kj, N, v);

					out[i][j].position.x += c_points[ki][kj].x * bi *bj;
					out[i][j].position.y += c_points[ki][kj].y * bi *bj;
					out[i][j].position.z += c_points[ki][kj].z * bi *bj;


					out[i][j].uv = glm::vec2(u,v);
				}
			}
		}
	}


	std::vector<Vertex> output;
	for (int i = 1; i < out.size(); i++)
	{
		const std::vector<Vertex>& prev = out[i-1];
		const std::vector<Vertex>& cur  = out[i-0];

		for (int j = 1; j < out[i].size(); j++) 
		{
			const Vertex& UL = cur[j - 1];
			const Vertex& UR = cur[j - 0];

			const Vertex& LR = prev[j - 0];
			const Vertex& LL = prev[j - 1];

			const glm::vec3 normal0 = glm::normalize(glm::cross(UR.position - LL.position, UL.position - LL.position));
			output.emplace_back(LL.position, normal0, LL.uv);
			output.emplace_back(UR.position, normal0, UR.uv);
			output.emplace_back(UL.position, normal0, UL.uv);

			const glm::vec3 normal1 = glm::normalize(glm::cross(LR.position - LL.position, UL.position - LL.position));
			output.emplace_back(LL.position, normal1, LL.uv);
			output.emplace_back(LR.position, normal1, LR.uv);
			output.emplace_back(UR.position, normal1, UR.uv);
		}
	}


	/*
	const std::vector<Vertex>& prev = out[out.size() - 1];
	const std::vector<Vertex>& cur = out[0];

	for (int j = 1; j < N_res; j++)
	{
		const Vertex& UL = cur[j - 1];
		const Vertex& UR = cur[j - 0];

		const Vertex& LR = prev[j - 0];
		const Vertex& LL = prev[j - 1];

		const glm::vec3 normal0 = glm::normalize(glm::cross(UR.position - LL.position, UL.position - LL.position));
		output.emplace_back(LL.position, normal0, LL.uv);
		output.emplace_back(UR.position, normal0, UR.uv);
		output.emplace_back(UL.position, normal0, UL.uv);

		const glm::vec3 normal1 = glm::normalize(glm::cross(LR.position - LL.position, UL.position - LL.position));
		output.emplace_back(LL.position, normal1, LL.uv);
		output.emplace_back(LR.position, normal1, LR.uv);
		output.emplace_back(UR.position, normal1, UR.uv);
	}*/
	

	//for (int i = 0; i < out.size(); i++)
	//{
	//	for (int j = 0; j < out[i].size(); j++)
	//	{
	//		output.emplace_back(out[i][j]);
	//	}
	//}


	for (int i = 0; i < output.size(); i++)
	{
		output[i].normal = (glm::dot(output[i].normal, output[i].position) < 0.0) ? output[i].normal : -output[i].normal;
	}


	return output;
}

//Carrega e linka shaders 
void loadShaders(void)
{
	AXIS_PROGRAM.compileShaders("axis.vertex" , "axis.frag");
	AXIS_PROGRAM.addAttribute("aPosition");
	AXIS_PROGRAM.linkShaders();

	POINTS_PROGRAM.compileShaders("point.vertex", "point.frag");
	POINTS_PROGRAM.addAttribute("pPosition");
	POINTS_PROGRAM.linkShaders();

	CURVE_PROGRAM.compileShaders("curve.vertex", "curve.frag");
	CURVE_PROGRAM.addAttribute("cPosition");
	CURVE_PROGRAM.linkShaders();

	AXIS3D_PROGRAM.compileShaders("axis3d.vertex", "axis3d.frag");
	AXIS3D_PROGRAM.addAttribute("a3Position");
	AXIS3D_PROGRAM.addAttribute("a3Color");
	AXIS3D_PROGRAM.linkShaders();

	SOLID_PROGRAM.compileShaders("solid.vertex", "solid2.frag");
	SOLID_PROGRAM.addAttribute("sPosition");
	SOLID_PROGRAM.addAttribute("sNormal");
	SOLID_PROGRAM.addAttribute("UV");
	SOLID_PROGRAM.linkShaders();

	

}

//VAO e VBO
void initVBOS(void)
{
	//Inaugurando uso dos vertex array objects (VAO), aumento de performance e finalmente entendi o que fazem de fato
	//VAOS controlam VBO's, permitem a reutilizacao de attribute locations em shaders diferentes, facilitam o agrupamento de objetos com caracteristicas em comum
	//Conferem ganho de velocidade de aprox 15%~ para casos genericos

	//axis
	glGenVertexArrays(1, &axisVAO);
	glBindVertexArray(axisVAO);
	glGenBuffers(1, &axisVBO);
	glBindBuffer(GL_ARRAY_BUFFER, axisVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(axis_data), axis_data, GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, axisVBO);
	glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,(void*)0);
	glDisableVertexAttribArray(0);
		

	//points
	glGenVertexArrays(1, &pointsVAO);
	glBindVertexArray(pointsVAO);
	glGenBuffers(1, &pointsVBO);
	
	glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
	
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2) * 20, NULL, GL_DYNAMIC_DRAW);
	
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
	glDisableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	//curve
	glGenVertexArrays(1, &curveVAO);
	glBindVertexArray(curveVAO);
	glGenBuffers(1, &curveVBO);
	glBindBuffer(GL_ARRAY_BUFFER, curveVBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2) * 2000, NULL, GL_DYNAMIC_DRAW);

	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, curveVBO);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
	glDisableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	
	//solid
	glGenVertexArrays(1, &solidVAO);
	glBindVertexArray(solidVAO);
	glGenBuffers(1, &solidVBO);
	glBindBuffer(GL_ARRAY_BUFFER, solidVBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * 30000, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, uv));

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	

	//axis3D
	glGenVertexArrays(1, &axis3DVAO);
	glBindVertexArray(axis3DVAO);
	glGenBuffers(1, &axis3DVBO);
	glBindBuffer(GL_ARRAY_BUFFER, axis3DVBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(axis_data2), axis_data2, GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, (0), (GLvoid*)0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, (0), (GLvoid*)(sizeof(float)*18));

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	
	
		
}

//Upload da textura
void genTex(void)
{
	texID = SOIL_load_OGL_texture("homog.png", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID , 0);


}

//calcula bounding box do modelo
void calcBB(const std::vector<Vertex>& _vert)
{

	//Obs muito mais simples do que eu imaginava, essencial para nao estragar a visao
	if (_vert.empty())
	{
		Bounds.center = glm::vec3(0.0, 0.0, 0.0);
		Bounds.min = glm::vec3(-300.0f, -300.0f, -300.0f);
		Bounds.max = glm::vec3(300.0f, 300.0f, 300.0f);
		return;
	}



	float min_x = _vert[0].position.x;
	float min_y = _vert[0].position.y;
	float min_z = _vert[0].position.z;

	float max_x = _vert[0].position.x;
	float max_y = _vert[0].position.y;
	float max_z = _vert[0].position.z;

	for (int i = 0; i < _vert.size(); i++)
	{
		//minimum
		if (_vert[i].position.x < min_x)
		{
			min_x = _vert[i].position.x;
		}


		if (_vert[i].position.y < min_y)
		{
			min_y = _vert[i].position.y;
		}


		if (_vert[i].position.z < min_z)
		{
			min_z = _vert[i].position.z;
		}


		//maximum
		if (_vert[i].position.x > max_x)
		{
			max_x = _vert[i].position.x;
		}


		if (_vert[i].position.y > max_y)
		{
			max_y = _vert[i].position.y;
		}


		if (_vert[i].position.z > max_z)
		{
			max_z = _vert[i].position.z;
		}

	}

	float centerX = max_x - ((max_x + fabs(min_x)) / 2);
	float centerY = max_y - ((max_y + fabs(min_y)) / 2);
	float centerZ = max_y - ((max_z + fabs(min_z)) / 2);

	Bounds.min = glm::vec3(min_x, min_y, min_z);
	Bounds.max = glm::vec3(max_x, max_y, max_z);
	Bounds.center = glm::vec3(centerX, centerY, centerZ);

}


/*-------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//MOUSE CALLBACK  muito importante
void onMouse(int button, int state, int x, int y)
{
	
	
		if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
		{

			if (x < WIDTH / 2)
			{
				control_points.emplace_back(glm::vec2(((float)x) - WIDTH / 4, HEIGHT / 2 - (float)y));

				if (is_parametric == false)
				{

					int size = control_points.size() - 1;
					p_curve.clear();

					float i = 0.0f;
					while (i <= 1.0f)
					{
						p_curve.emplace_back(evaluate_curve(control_points, i, size));
						i += 0.05f;
					}
					//std::cout << control_points.back().x << " " << control_points.back().y  << std::endl; debug

					p_curve.emplace_back(control_points.back()); //bresenham line

					model = genSOR(p_curve);
					sphericalCoords(model);
				}


				else if (is_parametric == true)
				{

					//rotatedCpoints.clear();
					rotatedCpoints.swap(rotate_points(control_points));
					//rotatedCpoints = rotate_points(control_points);
					//rotatedCpoints[rotatedCpoints.size() - 1] = rotatedCpoints[0];

					int size = control_points.size() - 1;
					p_curve.clear();

					float i = 0.0f;
					while (i <= 1.0f)
					{
						p_curve.emplace_back(evaluate_curve(control_points, i, size));
						i += 0.05f;
					}

					p_curve.emplace_back(control_points.back()); //bresenham line 

					model = bez_surface(rotatedCpoints);
				}
			}


			else if (x > WIDTH / 2)
			{
				arcball_on = true;
				last_mx = cur_mx = x;
				last_my = cur_my = y;
			}
		}
	

	



	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		
			if (x > WIDTH / 2)
			{
				arcball_on = false;
			}
				
	}



	if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
	{
		
			if (x < WIDTH / 2)
			{

				glm::vec2 point(x - (WIDTH/4), (HEIGHT/2) - y);
				for (int i = 0; i < control_points.size(); i++)
				{
					if (point.x >= (control_points[i].x - 5) && point.x <= (control_points[i].x + 5))
					{
						if (point.y >= (control_points[i].y - 5) && point.y <= (control_points[i].y + 5))
						{
							carry_point = true;
							c_index = i;

						}
					}
				}
			}
	}
		
	



	if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP)
	{
		
		if (x < WIDTH / 2)
		{
			carry_point = false;
			if (is_parametric == false)
			{
				model = genSOR(p_curve);
				sphericalCoords(model);
				calcBB(model); // remover caso performance caia
			}

			else if (is_parametric == true)
			{

				rotatedCpoints = rotate_points(control_points);
				model = bez_surface(rotatedCpoints);
				calcBB(model); // remover caso performance caia
				
			}
		}
	}


	
	glutPostRedisplay();
}

/*-------------------------------------------------------------------------------------------------------------------------------------------------------------*/

void onMotion(int x, int y)
{	
	if (arcball_on)
	{
		cur_mx = x;
		cur_my = y;

	}

	if(carry_point)
	{
		
		control_points[c_index] = glm::vec2(x - (WIDTH / 2)/2, HEIGHT / 2 - y);
		
		
		int size = control_points.size() - 1;
		p_curve.clear();

		float i = 0.0f;
		while (i <= 1.0f)
		{
			p_curve.emplace_back(evaluate_curve(control_points, i, size));
			i += 0.05f;
		}
		//std::cout << control_points.back().x << " " << control_points.back().y  << std::endl;

		p_curve.emplace_back(control_points.back()); //bresenham line

		//nota sobre otimização
		//model = genSOR(p_curve); extremamente pesado manter um novo solido a cada mouse_move mesmo em meu desktop (high end)

	}
		
	glutPostRedisplay();
}

//Calcula o vetor do eixo de rotacao
glm::vec3 get_arcball_vector(int _x, int _y)
{

	glm::vec3 P = glm::vec3((((_x - 600) / (WIDTH/4))) - 1.0, ((_y / (HEIGHT/2))) - 1.0, 0);
	
	//glm::vec3 P = glm::vec3(((_x / WIDTH)*2 ) - 1.0, ((_y / HEIGHT)*2)  - 1.0, 0);
	P.y = -P.y;
	float OP_squared = P.x * P.x + P.y * P.y;
	if (OP_squared <= 1)
	{
		P.z = sqrt(1 - OP_squared);
		P = glm::normalize(P);

	}
	else {
		P = glm::normalize(P);
	}
	return P;
}

void keyboard(uchar8 key, int x, int y)
{
	switch (key) 
	{
	case 27:
		exit(0);
	case 'c':
	case 'C':
		control_points.clear();
		p_curve.clear();
		model.clear();
		break;
	case '+':
		FOV -= 0.06f;
		break;
	case '-':
		FOV += 0.06f;
		break;
	case 'b':
	case 'B':
		calcBB(model);
		break;
	case 'p':
	case 'P':
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		break;
	case 'l':
	case 'L':
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		break;
	case 'F':
	case'f':
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case 'M':
	case 'm':
		
		break;
	case 'H':
	case 'h':
		
		break;
	case 'N':
	case 'n':
		
		break;
	case 'x':
	case 'X':
		see_X = true;
		break;
	case 'y':
	case 'Y':
		see_Y = true;
		break;
	case 'z':
	case 'Z':
		see_Z = true;
		break;
	case 'I':
	case 'i':
		see_Iso = true;
		break;
	case '1':
		texID = SOIL_load_OGL_texture("wood.png", SOIL_LOAD_AUTO, 1, 0);
		break;
	case '2' :
		texID = SOIL_load_OGL_texture("metal.png", SOIL_LOAD_AUTO, 1, 0);
		break;
	case '3':
		texID = SOIL_load_OGL_texture("homog.png", SOIL_LOAD_AUTO, 1, 0);
		break;
	case '0':

		fullscreen2D = false;
		if (fullscreen3D)
		{
			fullscreen3D = false;
		}
		else {

			fullscreen3D = true;
		}
		break;

	case '9':

		fullscreen3D = false;
		if (fullscreen2D)
		{
			fullscreen2D = false;
		}
		else {
			fullscreen2D = true;
		}
		break;
	case '8':

		fullscreen3D = false;
		fullscreen2D = false;
		break;
	case 'k':
	case 'K':
		is_parametric = true;
		break;
	case 'j':
	case 'J':
		is_parametric = false;
		break;
	default:
		printf("%c was pressed\n", key);
		break;


	}
	glutPostRedisplay();
}


void idleFunc(void)
{
	//
	

	glutPostRedisplay();
}


void reshape(int width, int height)
{
	//std::cout << "inside reshape" << std::endl;
	glViewport(0, 0, width, height);
	glutPostRedisplay();
}


void draw_axis(void)
{
	
		glViewport((GLsizei)0, (GLsizei)0, (GLsizei)WIDTH / 2, (GLsizei)HEIGHT);

				
		glBindVertexArray(axisVAO);

		AXIS_PROGRAM.use();

		glBindBuffer(GL_ARRAY_BUFFER, axisVBO);

		


		matrixID = AXIS_PROGRAM.getUniformLocation("matOrtho");
		glUniformMatrix4fv(matrixID, 1, GL_FALSE, &Ortho[0][0]);

		glDrawArrays(GL_LINES, 0, 4);

		AXIS_PROGRAM.unuse();
		//glBindVertexArray(0);
		//glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	
}


void draw_points(void)
{
	
		glViewport((GLsizei)0, (GLsizei)0, (GLsizei)WIDTH / 2, (GLsizei)HEIGHT);
				

		glBindVertexArray(pointsVAO);

		POINTS_PROGRAM.use();

		glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);

		


		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(glm::vec2)*control_points.size(), control_points.data());

		matrixID = POINTS_PROGRAM.getUniformLocation("pOrtho");
		glUniformMatrix4fv(matrixID, 1, GL_FALSE, &Ortho[0][0]);

		glDrawArrays(GL_POINTS, 0, control_points.size());

		POINTS_PROGRAM.unuse();
		//glBindVertexArray(0);
		//glBindBuffer(GL_ARRAY_BUFFER, 0);
	
}

//visualizacao da curva 2D (perfil)
void draw_curve(void)
{
	
		glViewport((GLsizei)0, (GLsizei)0, (GLsizei)WIDTH / 2, (GLsizei)HEIGHT);

		

		glBindVertexArray(curveVAO);

		CURVE_PROGRAM.use();

		glBindBuffer(GL_ARRAY_BUFFER, curveVBO);

		


		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(glm::vec2)*p_curve.size(), p_curve.data());

		matrixID = POINTS_PROGRAM.getUniformLocation("pOrtho");
		glUniformMatrix4fv(matrixID, 1, GL_FALSE, &Ortho[0][0]);


		glDrawArrays(GL_POINTS, 0, p_curve.size());

		POINTS_PROGRAM.unuse();

		//glBindVertexArray(0);
		//glBindBuffer(GL_ARRAY_BUFFER, 0);
	

}


void draw_axis3D(void)
{
	
	


	
		glViewport((GLsizei)WIDTH / 2, (GLsizei)0, (GLsizei)WIDTH / 2, (GLsizei)HEIGHT);


		glBindVertexArray(axis3DVAO);

		AXIS3D_PROGRAM.use();

		glBindBuffer(GL_ARRAY_BUFFER, axis3DVBO);



		matrixID = AXIS3D_PROGRAM.getUniformLocation("P");
		glUniformMatrix4fv(matrixID, 1, GL_FALSE, &Projection[0][0]);

		matrixID = AXIS3D_PROGRAM.getUniformLocation("V");
		glUniformMatrix4fv(matrixID, 1, GL_FALSE, &View[0][0]);

		matrixID = AXIS3D_PROGRAM.getUniformLocation("M");
		glUniformMatrix4fv(matrixID, 1, GL_FALSE, &Model[0][0]);

		glDrawArrays(GL_LINES, 0, 6);

		AXIS3D_PROGRAM.unuse();




}


void draw_solid(void)
{
	
		glViewport((GLsizei)WIDTH / 2, (GLsizei)0, (GLsizei)WIDTH / 2, (GLsizei)HEIGHT);


		if (!model.empty())
		{


			glBindVertexArray(solidVAO);

			SOLID_PROGRAM.use();

			glBindBuffer(GL_ARRAY_BUFFER, solidVBO);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, texID);





			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertex)*model.size(), model.data());

			//glm::mat4 MVP = Projection * View * Model;

			matrixID = SOLID_PROGRAM.getUniformLocation("P");
			glUniformMatrix4fv(matrixID, 1, GL_FALSE, &Projection[0][0]);

			matrixID = SOLID_PROGRAM.getUniformLocation("V");
			glUniformMatrix4fv(matrixID, 1, GL_FALSE, &View[0][0]);

			matrixID = SOLID_PROGRAM.getUniformLocation("M");
			glUniformMatrix4fv(matrixID, 1, GL_FALSE, &Model[0][0]);

			uniID = SOLID_PROGRAM.getUniformLocation("lightPos");
			glUniform3fv(uniID, 1, glm::value_ptr(lightPos));

			//uniID = SOLID_PROGRAM.getUniformLocation("lightColor");
			//glUniform3fv(uniID, 1, glm::value_ptr(lightColor));

			//uniID = SOLID_PROGRAM.getUniformLocation("objectColor");
			//glUniform3fv(uniID, 1 , glm::value_ptr(objectColor));

			uniID = SOLID_PROGRAM.getUniformLocation("eyePos");
			glUniform3fv(uniID, 1, glm::value_ptr(eyePos));

			//uniID = SOLID_PROGRAM.getUniformLocation("sStrength");
			//glUniform1f(uniID,  specular);

			uniID = SOLID_PROGRAM.getUniformLocation("material.ambient");
			glUniform3fv(uniID, 1, glm::value_ptr(current.ambient));
			uniID = SOLID_PROGRAM.getUniformLocation("material.diffuse");
			glUniform3fv(uniID, 1, glm::value_ptr(current.diffuse));
			uniID = SOLID_PROGRAM.getUniformLocation("material.specular");
			glUniform3fv(uniID, 1, glm::value_ptr(current.specular));
			uniID = SOLID_PROGRAM.getUniformLocation("material.shininess");
			glUniform1f(uniID, current.shininess);


			uniID = SOLID_PROGRAM.getUniformLocation("mysampler");
			glUniform1i(uniID, 0);

			glDrawArrays(GL_TRIANGLES, 0, model.size());

			SOLID_PROGRAM.unuse();

		}
	
}


void draw_text(void)
{
	if (!fullscreen3D && !fullscreen2D)
	{
		glWindowPos2i(50, 580);
		char text[] = "Curva de bezier 2D";
		glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text);


		glWindowPos2i(630, 565);
		glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)"X = ");

		glWindowPos2i(630, 550);
		glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)"Y = ");

		glWindowPos2i(630, 535);
		glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)"Z = ");




		if (is_parametric)
		{
			glWindowPos2i(650, 580);
			char text1[] = "SUPERFICIE DE BEZIER(PARAMETRICA)";
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text1);


			char minx[25];
			sprintf(minx, "%f", Bounds.min.x);

			char miny[25];
			sprintf(miny, "%f", Bounds.min.y);

			char minz[25];
			sprintf(minz, "%f", Bounds.min.z);

			glWindowPos2i(650, 565);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)minx);

			glWindowPos2i(650, 550);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)miny);

			glWindowPos2i(650, 535);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)minz);


			char cx[25];
			sprintf(cx, "%f", Bounds.center.x);

			char cy[25];
			sprintf(cy, "%f", Bounds.center.y);

			char cz[25];
			sprintf(cz, "%f", Bounds.center.z);

			glWindowPos2i(740, 565);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)cx);

			glWindowPos2i(740, 550);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)cy);

			glWindowPos2i(740, 535);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)cz);


			char maxx[25];
			sprintf(maxx, "%f", Bounds.max.x);

			char maxy[25];
			sprintf(maxy, "%f", Bounds.max.y);

			char maxz[25];
			sprintf(maxz, "%f", Bounds.max.z);

			glWindowPos2i(810, 565);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)maxx);

			glWindowPos2i(810, 550);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)maxy);

			glWindowPos2i(810, 535);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)maxz);






		}
		else {

			glWindowPos2i(650, 580);
			char text1[] = "SOLIDO DE REVOLUCAO ";
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text1);


			char minx[25];
			sprintf(minx, "%f", Bounds.min.x);

			char miny[25];
			sprintf(miny, "%f", Bounds.min.y);

			char minz[25];
			sprintf(minz, "%f", Bounds.min.z);

			glWindowPos2i(650, 565);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)minx);

			glWindowPos2i(650, 550);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)miny);

			glWindowPos2i(650, 535);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)minz);


			char cx[25];
			sprintf(cx, "%f", Bounds.center.x);

			char cy[25];
			sprintf(cy, "%f", Bounds.center.y);

			char cz[25];
			sprintf(cz, "%f", Bounds.center.z);

			glWindowPos2i(740, 565);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)cx);

			glWindowPos2i(740, 550);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)cy);

			glWindowPos2i(740, 535);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)cz);


			char maxx[25];
			sprintf(maxx, "%f", Bounds.max.x);

			char maxy[25];
			sprintf(maxy, "%f", Bounds.max.y);

			char maxz[25];
			sprintf(maxz, "%f", Bounds.max.z);

			glWindowPos2i(810, 565);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)maxx);

			glWindowPos2i(810, 550);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)maxy);

			glWindowPos2i(810, 535);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)maxz);


		}
	}
}


void draw(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//PLANO 2D
	draw_axis();
	
	//PONTOS
	draw_points();
	
	//CURVE
	draw_curve();

	float Max = std::max(Bounds.max.x, std::max(Bounds.max.y, Bounds.max.z));
	
	Projection = glm::perspective(FOV, (WIDTH / 2) / (HEIGHT), 0.01f, 20.0f*Bounds.max.z);

	//View = glm::lookAt(glm::vec3(1.5f*Max, 1.5f*Max, 1.5f*Max),
	//	glm::vec3(Bounds.center.x, Bounds.center.y, Bounds.center.z),
	//	glm::vec3(0.0f, 1.0f, 0.0f));

	
	if (cur_mx != last_mx || cur_my != last_my)
	{
		glm::vec3 va = get_arcball_vector(last_mx, last_my);
		glm::vec3 vb = get_arcball_vector(cur_mx, cur_my);

		float angle = acos(std::min(1.0f, glm::dot(va ,vb)));
				
		angle = angle*0.05;
		glm::vec3 axis_in_camera =  glm::cross(va, vb);
		
		
		glm::mat3 teste = glm::inverse(glm::mat3(View) * glm::mat3(Model));
		
		
		glm::vec3 axis_model = teste * axis_in_camera;

		

		//View = glm::rotate(View, glm::degrees(angle), axis_in_camera);

		//View = glm::rotate(View, glm::degrees(angle), axis_model);


		//Obs: Controle de camera é muito confuso

		//Model = glm::rotate(Model, glm::degrees(angle), axis_in_camera);

		Model = glm::rotate(Model, glm::degrees(angle), axis_model);
		
		last_mx = cur_mx;
		last_my = cur_my;

	}

	if (see_X)
	{
		View = glm::lookAt(glm::vec3(2.5f*Max, 0.0f*Max, 0.0f*Max),
				glm::vec3(Bounds.center.x, Bounds.center.y, Bounds.center.z),
				glm::vec3(0.0f, 1.0f, 0.0f));

		
		see_X = false;
	}
	
	if (see_Y)
	{
		View = glm::lookAt(glm::vec3(0.0f*Max, 2.5f*Max, 0.0f*Max),
			glm::vec3(Bounds.center.x, Bounds.center.y, Bounds.center.z),
			glm::vec3(1.0f, 0.0f, 0.0f));
		see_Y = false;

	}

	if (see_Z)
	{
		View = glm::lookAt(glm::vec3(0.0f*Max, 0.0f*Max, 2.5f*Max),
			glm::vec3(Bounds.center.x, Bounds.center.y, Bounds.center.z),
			glm::vec3(0.0f, 1.0f, 0.0f));
		see_Z = false;
	}

	if (see_Iso)
	{
		View = glm::lookAt(glm::vec3(2.5f*Max, 2.5f*Max, 2.5f*Max),
			glm::vec3(Bounds.center.x, Bounds.center.y, Bounds.center.z),
			glm::vec3(0.0f, 1.0f, 0.0f));
		see_Iso = false;
	}


	//Projection = glm::perspective(FOV, (WIDTH / 2) / (HEIGHT), 0.01f, 20.0f*Bounds.max.z);
	
	//View = glm::lookAt(glm::vec3(1.5f*Max, 1.5f*Max, 1.5f*Max),
	//				   glm::vec3(Bounds.center.x, Bounds.center.y, Bounds.center.z),
	//	               glm::vec3(0.0f, 1.0f, 0.0f));
	
	//axis3d
	draw_axis3D();
	
	//SOLID
	draw_solid();
	
	//Text
	draw_text();
	
	glutSwapBuffers();
}


///////////////////
///Inicialização///
///////////////////
void initGL(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize((int)WIDTH, (int)HEIGHT);
	glutCreateWindow("MATE71");

	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		std::cout << "could not load glew" << std::endl;
		std::cout << glewGetErrorString(err) << std::endl;
	}


	glutDisplayFunc(draw);
	glutReshapeFunc(reshape);

	glutKeyboardFunc(keyboard);
	glutIdleFunc(idleFunc);

	glutMouseFunc(onMouse);
	glutMotionFunc(onMotion);

	glEnable(GL_POINT_SPRITE);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
		
	

	//glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	//glCullFace(GL_BACK);
	glFrontFace(GL_CCW);

	glClearColor(0.0, 0.0, 0.0, 1.0);
	
}

int main(int argc, char **argv)
{
	initGL(&argc, argv);
	initVBOS();
	loadShaders();
	calcBB(model);

	genTex();

	float Max = std::max(Bounds.max.x, std::max(Bounds.max.y, Bounds.max.z));
	Projection = glm::perspective(FOV, (WIDTH / 2) / (HEIGHT), 0.01f, 20.0f*Bounds.max.z);

	View = glm::lookAt(glm::vec3(1.5f*Max, 1.5f*Max, 1.5f*Max),
		glm::vec3(Bounds.center.x, Bounds.center.y, Bounds.center.z),
		glm::vec3(0.0f, 1.0f, 0.0f));

	std::cout << "Press + or - to control FOV" << std::endl;
	std::cout << "Press c/C to clear profile curve" << std::endl;
	std::cout << "Click and Drag the mouse over the secondary window to move the camera" << std::endl;
	std::cout << "Press B to calculate the solid's bounding box " << std::endl;
	std::cout << "Press L to render solid in wireframe" << std::endl;
	std::cout << "Press P to render solid as points\n " << std::endl;

	std::cout << "Press X for a view on the X axis" << std::endl;
	std::cout << "Press Y for a view on the Y axis" << std::endl;
	std::cout << "Press Z for a view on the Z axis" << std::endl;
	std::cout << "Press I for an Isometric view\n" << std::endl;
	
	std::cout << "Press K for parametric surface" << std::endl;
	std::cout << "Press J for common surface" << std::endl;


	std::cout << "Press 1 for Wood Texture" << std::endl;
	std::cout << "Press 2 for Metal Texture" << std::endl;
	std::cout << "Press 3 for Homogeneous Texture" << std::endl;

	std::cout << "Press 0 for Fullscreen 3D" << std::endl;
	std::cout << "Press 9 for Fullscreen 2D" << std::endl;
	std::cout << "Press 8 to return" << std::endl;


	



	glutMainLoop();
	return 0;
}