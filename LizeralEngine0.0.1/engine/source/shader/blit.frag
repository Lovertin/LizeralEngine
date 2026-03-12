// shader/blit.frag
#version 460
layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D samplerColor;

void main() {
    // 纯粹的搬运，不加任何处理
    outColor = texture(samplerColor, inUV);
}