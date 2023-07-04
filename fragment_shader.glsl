/*
#version 330 core
out vec4 FragColor;
  
in vec4 vertexColor; // the input variable from the vertex shader (same name and same type)  

void main()
{
    FragColor = vertexColor;
} 
*/
#version 330 core

in vec2 TexCoord;
in vec2 LightIntensity;

out vec4 FragColor;

uniform sampler2D textureSampler;

void main()
{
    vec4 texColor = texture(textureSampler, TexCoord);
    float brighter = LightIntensity.x;
    if (LightIntensity.y > LightIntensity.x)
        brighter = LightIntensity.y;
    FragColor = texColor * vec4(vec3(brighter),1.0);
    //FragColor = vec4(1.0, 0,0,1.0);
}
