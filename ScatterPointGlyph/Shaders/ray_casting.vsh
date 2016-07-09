//[VERTEX SHADER]
#version 130

uniform mat4 projection, modelView;

in vec3 in_Position, in_Color;
out vec3 texCoord;

void main()
{
    texCoord = in_Position;
    gl_Position = projection * modelView * vec4(in_Position, 1.0);
}