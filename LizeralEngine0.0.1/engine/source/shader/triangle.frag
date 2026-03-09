#version 460

// 接收来自顶点着色器的颜色 (光栅化阶段会自动帮我们做平滑插值！)
layout(location = 0) in vec3 fragColor;

// 最终输出到 Swapchain 图片里的颜色
layout(location = 0) out vec4 outColor;

void main() {
    outColor = vec4(fragColor, 1.0); // 1.0 代表不透明 (Alpha)
}