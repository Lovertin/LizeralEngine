#version 460
layout(location = 0) in vec2 inUV;

layout(location = 0) out vec4 outGIHistory;    // History Length
layout(location = 1) out vec2 outMoments; 

layout(binding = 0) uniform sampler2D samplerNoisyGI;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;
layout(binding = 3) uniform sampler2D samplerVelocity;
layout(binding = 4) uniform sampler2D samplerPrevGIHistory;
layout(binding = 5) uniform sampler2D samplerPrevMoments;

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
    
    vec2 velocity = texture(samplerVelocity, inUV).xy;
    vec2 prevUV = inUV - velocity;

    bool isHistoryValid = (prevUV.x >= 0.0 && prevUV.x <= 1.0 && prevUV.y >= 0.0 && prevUV.y <= 1.0);
    // TODO: Future: add previous frame depth & normal validation to prevent ghosting from occluded objects. Currently using UV validation

    float historyLength = 0.0;
    vec3 prevGI = vec3(0.0);
    vec2 prevMoments = vec2(0.0);

    if (isHistoryValid) {
        vec4 prevGIHistory = texture(samplerPrevGIHistory, prevUV);
        prevGI = prevGIHistory.rgb;
        historyLength = prevGIHistory.a;
        prevMoments = texture(samplerPrevMoments, prevUV).rg;
    }

    historyLength = min(64.0, isHistoryValid ? historyLength + 1.0 : 1.0);

    // Calculate Blending Weight (alpha)
    // For newly appearing objects: alpha = 1 (only trust current frame) For stable objects: alpha = 0.05 (trust history more, extremely smooth)
    float alpha = max(0.02, 1.0 / historyLength);

    // Temporal Blending GI Color
    vec3 resolvedGI = mix(prevGI, currentGI, alpha);
    
    // Temporal Blending Variance (Moments)
    vec2 currentMoments = vec2(currentLuma, currentLuma * currentLuma);
    vec2 resolvedMoments = mix(prevMoments, currentMoments, alpha);

    outGIHistory = vec4(resolvedGI, historyLength);
    outMoments = resolvedMoments;
}