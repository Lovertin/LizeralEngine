#version 460
layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D samplerColor;

// 把 ACESFilm 搬到这里来！
vec3 ACESFilm(vec3 x) {
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0, 1.0);
}

void main() {
    // 1. 读取线性空间的纯净 HDR 画面 (这是 TAA 处理后的结果)
    vec3 hdrColor = texture(samplerColor, inUV).rgb;
    
    float exposure = 1.0; 
    hdrColor *= exposure;
    
    // 2. 进行 ToneMapping 映射到屏幕支持的 LDR 色彩
    vec3 mappedColor = ACESFilm(hdrColor);
    
    // 3. Gamma 校正 (极其重要！屏幕是以 sRGB 显示的，我们需要转一下)
    // 很多时候你觉得画面暗或者对比度怪，就是因为漏了 Gamma 2.2 校正
    mappedColor = pow(mappedColor, vec3(1.0 / 2.2));

    outColor = vec4(mappedColor, 1.0);
}