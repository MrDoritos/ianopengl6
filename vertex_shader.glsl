/*#version 330 core
layout (location = 0) in vec3 aPos; // the position variable has attribute position 0
  
out vec4 vertexColor; // specify a color output to the fragment shader

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

void main()
{
    gl_Position = proj * view * model * vec4(aPos, 1.0); // see how we directly give a vec3 to vec4's constructor
    vertexColor = vec4(0.5, 0.0, 0.0, 1.0); // set the output variable to a dark-red color
}

*/
#version 330 core

layout (location = 0) in vec3 aPosition;
layout (location = 1) in vec2 aTexCoord;
layout (location = 2) in vec2 aLightIntensity;

out vec2 TexCoord;
out vec2 LightIntensity;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

void main()
{
    
    gl_Position = proj * view * model * vec4(aPosition, 1.0);
    TexCoord = aTexCoord;
    LightIntensity = aLightIntensity;
}
