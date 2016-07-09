#version 330

uniform mat4 ModelViewProj;
uniform vec3 VertexScale;

in vec4 PositionIn;
out vec4 Position;

void main() 
{
	gl_Position = ModelViewProj * vec4(PositionIn.xyz * VertexScale - VertexScale * 0.5, 1.0);
    Position = PositionIn;
}