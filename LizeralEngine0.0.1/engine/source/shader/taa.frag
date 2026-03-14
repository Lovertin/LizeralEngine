#version 460

layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor; 

// ★ 贴图减少到了 3 张！
layout(binding = 0) uniform sampler2D samplerSceneColor; // Denoise传来的画面
layout(binding = 1) uniform sampler2D samplerHistory;
layout(binding = 2) uniform sampler2D samplerVelocity;

vec3 RGBToYCoCg(vec3 rgb) {
    return vec3(
         0.25 * rgb.r + 0.5 * rgb.g + 0.25 * rgb.b,
         0.5  * rgb.r               - 0.5  * rgb.b,
        -0.25 * rgb.r + 0.5 * rgb.g - 0.25 * rgb.b
    );
}

vec3 YCoCgToRGB(vec3 ycocg) {
    return vec3(
        ycocg.x + ycocg.y - ycocg.z,
        ycocg.x + ycocg.z,
        ycocg.x - ycocg.y - ycocg.z
    );
}

void main() {
    // 1. 直接读取当前帧最终画面
    vec4 currentColor = texture(samplerSceneColor, inUV);
    
    vec2 texelSize = 1.0 / textureSize(samplerSceneColor, 0);

    vec2 velocity = vec2(0.0);
    float maxVelocitySq = -1.0;
    for (int y = -1; y <= 1; ++y) {
        for (int x = -1; x <= 1; ++x) {
            vec2 v = texture(samplerVelocity, inUV + vec2(x, y) * texelSize).xy;
            float vSq = dot(v, v);
            if (vSq > maxVelocitySq) {
                maxVelocitySq = vSq;
                velocity = v;
            }
        }
    }
    
    vec2 prevUV = inUV - velocity;
    if (prevUV.x < 0.0 || prevUV.x > 1.0 || prevUV.y < 0.0 || prevUV.y > 1.0) {
        outColor = currentColor;
        return;
    }

    vec4 historyColor = texture(samplerHistory, prevUV);

    // 2. 方差裁剪 (现在只需要采样 1 张图了！性能起飞)
    vec3 m1 = vec3(0.0);
    vec3 m2 = vec3(0.0);

    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            // ★ 只采样一次 SceneColor！
            vec3 nColor = texture(samplerSceneColor, inUV + vec2(x, y) * texelSize).rgb;
            
            vec3 ycocg = RGBToYCoCg(nColor);
            m1 += ycocg;
            m2 += ycocg * ycocg;
        }
    }

    m1 /= 9.0;
    m2 /= 9.0;

    vec3 mu = m1;
    vec3 sigma = sqrt(max(m2 - m1 * m1, vec3(0.0)));

    float gamma = 1.2; // 裁剪强度
    vec3 colorMin = mu - gamma * sigma;
    vec3 colorMax = mu + gamma * sigma;

    vec3 historyYCoCg = RGBToYCoCg(historyColor.rgb);
    historyYCoCg = clamp(historyYCoCg, colorMin, colorMax);
    historyColor.rgb = YCoCgToRGB(historyYCoCg);

    float blendWeight = 0.05; 
    float velocityLength = length(velocity);
    if (velocityLength > 0.001) {
        blendWeight = mix(0.05, 0.3, clamp(velocityLength * 100.0, 0.0, 1.0));
    }

    outColor = mix(historyColor, currentColor, blendWeight);
}