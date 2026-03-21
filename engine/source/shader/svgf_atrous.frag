#version 460
//extension for BDA 
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require

layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outDenoisedGI;

layout(binding = 0) uniform sampler2D samplerGIHistory;
layout(binding = 1) uniform sampler2D samplerMoments;
layout(binding = 2) uniform sampler2D samplerNormalRoughness;
layout(binding = 3) uniform sampler2D samplerDepth;

// be the same with lighting.frag

struct GlobalFrameData {
    mat4 viewProj;
    mat4 invViewProj;
    mat4 prevViewProj;
    vec3 cameraPos;
    float padding1;
    vec3 lightDir;
    float lightIntensity;
    vec3 lightColor;
    uint frameIndex;
    vec2 jitter;
    vec2 padding2;
};

layout(buffer_reference, scalar, buffer_reference_align = 16) readonly buffer FrameDataBuffer {
    GlobalFrameData frame;
};

layout(push_constant) uniform PushConstants {
    FrameDataBuffer frameData;
    uint64_t instanceDescAddr;       
    uint64_t pointLightsAddr;
    uint stepSize;
    uint pointLightCount;
} pc;

vec3 ReconstructWorldPos(vec2 uv, float depth) {
    vec4 ndc = vec4(uv * 2.0 - 1.0, depth, 1.0);
    vec4 worldPos = pc.frameData.frame.invViewProj * ndc; 
    return worldPos.xyz / worldPos.w;
}

float GetLuminance(vec3 color) {
    return dot(color, vec3(0.2126, 0.7152, 0.0722));
}

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

    float variance = max(0.0, moments.y - moments.x * moments.x);
    
    int stepSize = int(pc.stepSize); 

    float phiLuma = 10.0; 
    float phiNormal = 16.0;

    vec3 sumGI = centerGI;
    float sumWeight = 1.0;

    vec2 texSize = textureSize(samplerGIHistory, 0);
    vec2 texelSize = 1.0 / texSize;

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

            float wSpatial = kernelWeights[abs(x)] * kernelWeights[abs(y)];
            float wNormal = pow(max(dot(centerNormal, sampleNormal), 0.0), phiNormal);
            float wLuma = exp(-abs(centerLuma - sampleLuma) / (phiLuma * sqrt(variance) + 0.001));

            float weight = wSpatial * wNormal * wLuma;

            sumGI += sampleGI * weight;
            sumWeight += weight;
        }
    }

    outDenoisedGI = vec4(sumGI / sumWeight, 1.0);
}