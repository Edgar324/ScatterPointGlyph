#version 330

in vec4 PositionIn;
out vec3 Position;

void main() 
{
	gl_Position = vec4(PositionIn.xyz, 1.0);
    Position = vec3((PositionIn.xy + vec2(1.0, 1.0)) / 2.0, Position.z);
}