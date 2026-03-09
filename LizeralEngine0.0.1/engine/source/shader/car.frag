#version 460
layout(location = 0) in vec3 fragNormal;
layout(location = 0) out vec4 outColor;

void main() {
    // ★ 强制把法线转成 RGB 颜色，如果画出来了，绝对是色彩斑斓的，不可能黑屏！
    vec3 color = normalize(fragNormal) * 0.5 + 0.5; 
    outColor = vec4(color, 1.0);
}