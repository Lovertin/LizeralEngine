#version 460
#extension GL_EXT_nonuniform_qualifier : require 
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require

layout(location = 0) in vec3 fragNormal;
layout(location = 1) in vec2 fragUV;
layout(location = 2) flat in uint fragMatID;
layout(location = 3) in vec4 fragCurrentPosClip;
layout(location = 4) in vec4 fragPrevPosClip;
layout(location = 5) in vec3 fragWorldPos;  

layout(location = 2) out vec2 outVelocity;
layout(location = 0) out vec4 outAlbedoMetallic; // RT0
layout(location = 1) out vec4 outNormalRoughness; // RT1

layout(binding = 0) uniform sampler2D GlobalTextures[1024]; 

struct Material {
    vec4 baseColorFactor;
    float metallicFactor;
    float roughnessFactor;
    int albedoTex;
    int normalTex;
    int ormTex;
    int emissiveTex;
    int pad0, pad1;
};

layout(buffer_reference, scalar, buffer_reference_align = 4) readonly buffer MaterialBuffer { Material m[]; };

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

layout(buffer_reference, scalar, buffer_reference_align = 16) readonly buffer InstanceDataBuffer {
    mat4 model;
    mat4 prevModel;
};

layout(push_constant) uniform PushConstants {
    FrameDataBuffer frameDataAddr;       // 8 words
    InstanceDataBuffer instanceDataAddr; // 8 words
    
    uint64_t vBuf;
    uint64_t mBuf;
    uint64_t iBuf;
    uint64_t bBuf;
    uint64_t matBuf;
    uint totalMeshlets;
    uint textureOffset;
} pc;

vec3 PerturbNormal(vec3 worldNormal, vec3 worldPos, vec2 uv, vec3 normalMapSample) {
    vec3 q1  = dFdx(worldPos);
    vec3 q2  = dFdy(worldPos);
    vec2 st1 = dFdx(uv);
    vec2 st2 = dFdy(uv);

    vec3 N   = normalize(worldNormal);
    vec3 T   = normalize(q1 * st2.t - q2 * st1.t);
    vec3 B   = -normalize(cross(N, T));
    mat3 TBN = mat3(T, B, N);

    vec3 tNormal = normalMapSample * 2.0 - 1.0;
    return normalize(TBN * tNormal);
}

void main() {
    MaterialBuffer matBuf = MaterialBuffer(pc.matBuf);
    uint localMatID = fragMatID - pc.textureOffset; 
    Material mat = matBuf.m[localMatID];

    vec4 albedo = mat.baseColorFactor;
    if (mat.albedoTex >= 0) {
        albedo *= texture(GlobalTextures[nonuniformEXT(mat.albedoTex)], fragUV);
    }

    float ao = 1.0;
    float roughness = mat.roughnessFactor;
    float metallic = mat.metallicFactor;
    
    if (mat.ormTex >= 0) {
        vec3 orm = texture(GlobalTextures[nonuniformEXT(mat.ormTex)], fragUV).rgb;
        ao = orm.r;              // R :AO
        roughness *= orm.g;      // G :Roughness
        metallic *= orm.b;       // B :Metallic
    }

    vec3 finalNormal = normalize(fragNormal);
    if (mat.normalTex >= 0) {
        vec3 normalMapSample = texture(GlobalTextures[nonuniformEXT(mat.normalTex)], fragUV).rgb;
        finalNormal = PerturbNormal(finalNormal, fragWorldPos, fragUV, normalMapSample);
    }

    vec2 currentNDC = fragCurrentPosClip.xy / fragCurrentPosClip.w;
    vec2 prevNDC = fragPrevPosClip.xy / fragPrevPosClip.w;
    vec2 currentUV = currentNDC * 0.5 + 0.5;
    vec2 prevUV = prevNDC * 0.5 + 0.5;
    outVelocity = currentUV - prevUV;

    vec3 finalAlbedo = albedo.rgb * ao;

    // RT0: RGB save color as Alpha save metallic
    outAlbedoMetallic = vec4(finalAlbedo, metallic);

    // RT1: RGB save normal ,Alpha save roughness
    outNormalRoughness = vec4(finalNormal, roughness);
}