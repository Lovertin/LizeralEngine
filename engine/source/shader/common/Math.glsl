// shader/Common/Math.glsl
#ifndef MATH_GLSL
#define MATH_GLSL

#define PI 3.14159265359

uint pcg_hash(uint seed) {
    uint state = seed * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

float rand(inout uint seed) {
    seed = pcg_hash(seed);
    return float(seed & 0x00FFFFFFu) / 16777216.0;
}

// 解耦：通过参数传入 invViewProj，而不是硬编码读取全局 pc
vec3 ReconstructWorldPos(vec2 uv, float depth, mat4 invViewProj) {
    vec4 ndc = vec4(uv * 2.0 - 1.0, depth, 1.0);
    vec4 worldPos = invViewProj * ndc; 
    return worldPos.xyz / worldPos.w;
}

// UE4 平方反比衰减算法 (也从主函数中剥离出来)
float CalculatePointLightFalloff(float distance, float radius) {
    float d_over_r = distance / radius;
    float d_over_r_4 = d_over_r * d_over_r * d_over_r * d_over_r;
    float falloff = clamp(1.0 - d_over_r_4, 0.0, 1.0);
    return (falloff * falloff) / (distance * distance + 1.0);
}

#endif