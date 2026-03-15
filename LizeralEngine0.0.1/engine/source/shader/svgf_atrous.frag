#version 460
layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outDenoisedGI;

layout(binding = 0) uniform sampler2D samplerGIHistory;
layout(binding = 1) uniform sampler2D samplerMoments;
layout(binding = 2) uniform sampler2D samplerNormalRoughness;
layout(binding = 3) uniform sampler2D samplerDepth;

layout(push_constant) uniform PushConstants {
    mat4 invViewProj;
    mat4 viewProj;
    vec3 cameraPos;
    float padding;
    uint frameIndex;
    uint stepSize;
} pc;

vec3 ReconstructWorldPos(vec2 uv, float depth) {
    vec4 ndc = vec4(uv * 2.0 - 1.0, depth, 1.0);
    vec4 worldPos = pc.invViewProj * ndc; 
    return worldPos.xyz / worldPos.w;
}

float GetLuminance(vec3 color) {
    return dot(color, vec3(0.2126, 0.7152, 0.0722));
}

// 3x3 离散高斯模糊核权重
const float kernelWeights[3] = float[](1.0, 2.0/3.0, 1.0/6.0);

void main() {
    float centerDepth = texture(samplerDepth, inUV).r;
    if (centerDepth <= 0.0) {
        outDenoisedGI = texture(samplerGIHistory, inUV);
        return;
    }

    vec3 centerGI = texture(samplerGIHistory, inUV).rgb;
    vec2 moments = texture(samplerMoments, inUV).rg;
    vec3 centerNormal = normalize(texture(samplerNormalRoughness, inUV).rgb + vec3(0.0001));
    vec3 centerWorldPos = ReconstructWorldPos(inUV, centerDepth);
    float centerLuma = GetLuminance(centerGI);

    // 1. 核心公式：计算方差 Variance = E(x^2) - E(x)^2
    float variance = max(0.0, moments.y - moments.x * moments.x);
    
    // 我们在这个 Pass 写死步长为 1。如果要实现极致的 SVGF，
    // 我们应该在 C++ 里 Dispatch 这个 Shader 3 次，每次传入 stepSize = 1, 2, 4
    int stepSize = int(pc.stepSize); 

    // 控制各项权重的灵敏度参数
    float phiLuma = 4.0;
    float phiNormal = 128.0;

    vec3 sumGI = centerGI;
    float sumWeight = 1.0;

    vec2 texSize = textureSize(samplerGIHistory, 0);
    vec2 texelSize = 1.0 / texSize;

    // 5x5 À-Trous 迭代
    for (int y = -2; y <= 2; ++y) {
        for (int x = -2; x <= 2; ++x) {
            if (x == 0 && y == 0) continue;

            vec2 offset = vec2(x, y) * float(stepSize);
            vec2 sampleUV = inUV + offset * texelSize;

            if (sampleUV.x < 0.0 || sampleUV.x > 1.0 || sampleUV.y < 0.0 || sampleUV.y > 1.0) continue;

            float sampleDepth = texture(samplerDepth, sampleUV).r;
            if (sampleDepth <= 0.0) continue;

            vec3 sampleGI = texture(samplerGIHistory, sampleUV).rgb;
            vec3 sampleNormal = normalize(texture(samplerNormalRoughness, sampleUV).rgb + vec3(0.0001));
            float sampleLuma = GetLuminance(sampleGI);

            // 权重 1: 空间高斯权重
            float wSpatial = kernelWeights[abs(x)] * kernelWeights[abs(y)];

            // 权重 2: 法线相似度权重 (越接近 1 说明越平坦)
            float wNormal = pow(max(dot(centerNormal, sampleNormal), 0.0), phiNormal);

            // 权重 3: 亮度差异权重 (★ 方差引导的核心！)
            // 方差(variance)越大，分母越大，wLuma 越接近 1 (狠狠地模糊)
            // 方差(variance)越小，分母越小，哪怕亮度差一点点，wLuma 也会骤降 (保留锐利边缘)
            float wLuma = exp(-abs(centerLuma - sampleLuma) / (phiLuma * sqrt(variance) + 0.001));

            float weight = wSpatial * wNormal * wLuma;

            sumGI += sampleGI * weight;
            sumWeight += weight;
        }
    }

    outDenoisedGI = vec4(sumGI / sumWeight, 1.0);
}