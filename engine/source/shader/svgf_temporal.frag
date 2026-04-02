#version 460
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require
layout(location = 0) in vec2 inUV;

layout(location = 0) out vec4 outGIHistory;    // History Length
layout(location = 1) out vec2 outMoments; 

layout(binding = 0) uniform sampler2D samplerNoisyGI;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;
layout(binding = 3) uniform sampler2D samplerVelocity;
layout(binding = 4) uniform sampler2D samplerPrevGIHistory;
layout(binding = 5) uniform sampler2D samplerPrevMoments;

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
    FrameDataBuffer frameDataAddr;
    uint64_t instanceDescAddr;
    uint64_t pointLightsAddr;
    uint stepSize;
    uint pointLightCount;
} pc;

float GetLuminance(vec3 color) {
    return dot(color, vec3(0.2126, 0.7152, 0.0722));
}

void main() {
    float depth = texture(samplerDepth, inUV).r;
    if (depth <= 0.0) {
        outGIHistory = vec4(texture(samplerNoisyGI, inUV).rgb, 0.0);
        outMoments = vec2(0.0);
        return;
    }

    vec3 currentGI = texture(samplerNoisyGI, inUV).rgb;
    float currentLuma = GetLuminance(currentGI);
    if (pc.frameDataAddr.frame.frameIndex < 3u) {
        outGIHistory = vec4(currentGI, 1.0);
        outMoments = vec2(currentLuma, currentLuma * currentLuma);
        return;
    }
    
    vec2 velocity = texture(samplerVelocity, inUV).xy;
    vec2 prevUV = inUV - velocity;
    float velocityLen = length(velocity);

    bool isHistoryValid = (prevUV.x >= 0.0 && prevUV.x <= 1.0 && prevUV.y >= 0.0 && prevUV.y <= 1.0);
    // Reject history for unstable reprojection to avoid fullscreen flicker after mode/resolution changes.
    isHistoryValid = isHistoryValid && (velocityLen < 0.5);

    float historyLength = 0.0;
    vec3 prevGI = vec3(0.0);
    vec2 prevMoments = vec2(0.0);

    if (isHistoryValid) {
        vec4 prevGIHistory = texture(samplerPrevGIHistory, prevUV);
        float prevHistoryLength = clamp(prevGIHistory.a, 0.0, 64.0);
        bool historyFinite = !any(isnan(prevGIHistory.rgb)) && !any(isinf(prevGIHistory.rgb));
        isHistoryValid = (prevHistoryLength > 0.5) && historyFinite;
        if (isHistoryValid) {
            prevGI = prevGIHistory.rgb;
            historyLength = prevHistoryLength;
            prevMoments = texture(samplerPrevMoments, prevUV).rg;
        }
    }

    historyLength = min(64.0, isHistoryValid ? historyLength + 1.0 : 1.0);

    // Calculate Blending Weight (alpha)
    // For newly appearing objects: alpha = 1 (only trust current frame) For stable objects: alpha = 0.05 (trust history more, extremely smooth)
    float alpha = max(0.06, 1.0 / historyLength);
    alpha = max(alpha, clamp(velocityLen * 2.5, 0.0, 0.35));

    // Temporal Blending GI Color
    vec3 resolvedGI = mix(prevGI, currentGI, alpha);
    
    // Temporal Blending Variance (Moments)
    vec2 currentMoments = vec2(currentLuma, currentLuma * currentLuma);
    vec2 resolvedMoments = mix(prevMoments, currentMoments, alpha);

    outGIHistory = vec4(resolvedGI, historyLength);
    outMoments = resolvedMoments;
}
