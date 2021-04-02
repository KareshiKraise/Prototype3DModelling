#version 330

in vec3 o_a3Color;

void main()
{
	gl_FragColor = vec4(o_a3Color,1.0);
}