#version 330 core
layout (location = 0) in vec3 aPos;

out vec3 TexCoords; // 向片段着色器传递 3D 纹理采样坐标

uniform mat4 projection;
uniform mat4 view;

void main()
{
    TexCoords = aPos;
    
    vec4 pos = projection * view * vec4(aPos, 1.0);

    gl_Position = pos.xyww;
}