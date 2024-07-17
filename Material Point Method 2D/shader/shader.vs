#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 3) in mat4 aInstanceMatrix;

out vec3 fColor;

uniform mat4 projection;
uniform mat4 view;

void main()
{
    //fColor = vec3(0.0, 0.5, 0.5);
    fColor = vec3(0.8, 0.5, 0.0);
    //gl_Position = vec4(aPos + aOffset, 1.0);
    gl_Position = projection * view * aInstanceMatrix * vec4(aPos, 1.0f); 
}