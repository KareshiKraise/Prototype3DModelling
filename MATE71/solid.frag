#version 330


in vec3 FragPos;
in vec3 slNormal;

uniform vec3 lightPos;
uniform vec3 objectColor;
uniform vec3 lightColor;
uniform vec3 eyePos;

uniform float sStrength;

void main()
{	

	
	
	float ambientStrength = 0.25f;
	vec3 ambient = ambientStrength * objectColor;
	

	vec3 norm = normalize(slNormal);
	vec3 lightDir = normalize(lightPos - FragPos);


	float diff = max(dot(lightDir, norm), 0.0);
	vec3 diffuse = diff * objectColor;


	float specStrength = sStrength;
	vec3 viewDir = normalize(eyePos - FragPos);
	
	
	
	vec3 half = normalize(lightDir + viewDir);
	float spec = pow(max(dot(norm, half), 0.0), 32.0);

	vec3 specular = specStrength * lightColor * spec ;

	vec4 color = vec4(ambient + diffuse + specular, 1.0f);




	

	//vec3 specular = specStrength * spec * lightColor;
	//vec3 result = (ambient + diffuse + specular) * objectColor;
	//vec4 color = vec4(result, 1.0f);

	//first value is ambient
	//vec3 result = (vec3(0.15,0.15,0.15) + diffuse) * objectColor;
	//vec4 color = vec4(result, 1.0f);

	gl_FragColor = color;

	


}