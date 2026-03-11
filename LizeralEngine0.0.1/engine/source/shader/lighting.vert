#version 460

layout(location = 0) out vec2 outUV;

void main() {
    // 利用位运算生成全屏三角形的 UV 坐标 (0,0), (2,0), (0,2)
    outUV = vec2((gl_VertexIndex << 1) & 2, gl_VertexIndex & 2);
    // 映射到 Vulkan 的 NDC 坐标系 [-1, 1]
    gl_Position = vec4(outUV * 2.0f - 1.0f, 0.0f, 1.0f);
}