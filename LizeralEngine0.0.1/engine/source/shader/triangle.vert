#version 460

// 输出给下一个阶段（片元着色器）的颜色
layout(location = 0) out vec3 fragColor;

// 直接在 GPU 里硬编码 3 个顶点的位置 (屏幕坐标范围是 -1 到 1)
vec2 positions[3] = vec2[](
    vec2(0.0, -0.5),  // 上中
    vec2(0.5, 0.5),   // 右下
    vec2(-0.5, 0.5)   // 左下
);

// 硬编码 3 个顶点的颜色 (RGB)
vec3 colors[3] = vec3[](
    vec3(1.0, 0.0, 0.0), // 红
    vec3(0.0, 1.0, 0.0), // 绿
    vec3(0.0, 0.0, 1.0)  // 蓝
);

void main() {
    // gl_VertexIndex 是 Vulkan 自动传进来的编号 (0, 1, 2)
    gl_Position = vec4(positions[gl_VertexIndex], 0.0, 1.0);
    fragColor = colors[gl_VertexIndex];
}