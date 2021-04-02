#version 330

struct Material {
	vec3 ambient;
	vec3 diffuse;
	vec3 specular;
	float shininess;
};

in vec3 FragPos;
in vec3 slNormal;
in vec2 sUV;

uniform vec3 lightPos;
//uniform vec3 objectColor;
//uniform vec3 lightColor;
uniform vec3 eyePos;
uniform sampler2D mysampler;

uniform Material material;


void main()
{	

	vec3 color = texture(mysampler, sUV).rgb;

	//ambient
	vec3 ambient = color * material.ambient;

	//diffuse
	vec3 lightDir = normalize(lightPos - FragPos);
	vec3 normal = normalize(slNormal);
	float diff = max(dot(lightDir, normal), 0.0);
	vec3 diffuse = diff * material.diffuse * color;


	//specular
	vec3 viewDir = normalize(eyePos - FragPos);
	vec3 halfwayDir = normalize(lightDir + viewDir);
	float spec = pow(max(dot(normal, halfwayDir), 0.0), material.shininess);

	vec3 specular = spec * material.specular ;

	gl_FragColor = vec4(ambient + diffuse + specular, 1.0f);




	/*vec3 text = vec3(texture2D(mysampler, sUV));
	
	vec3 ambient = (material.ambient) * text;
	

	vec3 norm = normalize(slNormal);
	vec3 lightDir = normalize(lightPos - FragPos);


	float diff = max(dot(norm, lightDir), 0.0);
	vec3 diffuse = (diff * material.diffuse) * text;
	
	
	vec3 viewDir = normalize(eyePos - FragPos);
	
	
	vec3 half = normalize(lightDir + viewDir);
	float spec = pow(max(dot(norm, half), 0.0), material.shininess);

	vec3 specular = (spec * material.specular);

	vec4 color = vec4(ambient + diffuse + specular, 1.0f);
			

	gl_FragColor = color;*/
	//gl_FragColor = vec4(slNormal, 1.0);

	


}