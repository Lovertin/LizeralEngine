#version 460
#extension GL_EXT_nonuniform_qualifier : require 
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require

layout(location = 0) in vec3 fragNormal;
layout(location = 1) in vec2 fragUV;
layout(location = 2) flat in uint fragTexID;

layout(location = 3) in vec4 fragCurrentPosClip;
layout(location = 4) in vec4 fragPrevPosClip;

layout(location = 2) out vec2 outVelocity;

layout(location = 0) out vec4 outAlbedoMetallic; // RT0
layout(location = 1) out vec4 outNormalRoughness; // RT1

layout(binding = 0) uniform sampler2D GlobalTextures[1024]; 

struct Material {
    vec4 baseColorFactor;
    float metallicFactor;
    float roughnessFactor;
    vec2 padding;
};

layout(buffer_reference, scalar, buffer_reference_align = 4) readonly buffer MaterialBuffer { Material m[]; };

// ==========================================================
// ★ 1. 将旧的 pc 结构替换为全新的 BDA 指针架构
// ==========================================================
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

// 映射全局帧数据指针
layout(buffer_reference, scalar, buffer_reference_align = 16) readonly buffer FrameDataBuffer {
    GlobalFrameData frame;
};

// 映射实例矩阵数据指针
layout(buffer_reference, scalar, buffer_reference_align = 16) readonly buffer InstanceDataBuffer {
    mat4 model;
    mat4 prevModel;
};

// ★ 2. 瘦身版 PushConstants (必须和 car.mesh / car.task 100% 相同)
layout(push_constant) uniform PushConstants {
    FrameDataBuffer frameDataAddr;       // 8 字节：指向全局数据
    InstanceDataBuffer instanceDataAddr; // 8 字节：指向本物体的矩阵
    
    uint64_t vBuf;
    uint64_t mBuf;
    uint64_t iBuf;
    uint64_t bBuf;
    uint64_t matBuf;
    uint totalMeshlets;
    uint textureOffset;
} pc;

void main() {
    
    // 1. 获取材质数据
    MaterialBuffer matBuf = MaterialBuffer(pc.matBuf);
    // ★ 修复：减去当前模型的偏移量，拿到局部材质 ID
    uint localMatID = fragTexID - pc.textureOffset; 
    Material mat = matBuf.m[localMatID];

    // 2. 无绑定贴图采样
    uint texIndex = fragTexID % 1024; 
    vec4 texColor = texture(GlobalTextures[nonuniformEXT(texIndex)], fragUV);

    vec2 currentNDC = fragCurrentPosClip.xy / fragCurrentPosClip.w;
    vec2 prevNDC = fragPrevPosClip.xy / fragPrevPosClip.w;

    vec2 currentUV = currentNDC * 0.5 + 0.5;
    vec2 prevUV = prevNDC * 0.5 + 0.5;

    outVelocity = currentUV - prevUV;

    // 3. 计算 PBR 基础参数
    vec4 albedo = texColor * mat.baseColorFactor;

    // ★ 4. 暴力填入 G-Buffer，光照计算被彻底剥离到下一个 Pass！
    // RT0: RGB 存颜色，Alpha 存金属度
    outAlbedoMetallic = vec4(albedo.rgb, mat.metallicFactor);

    // RT1: RGB 存法线(记得归一化)，Alpha 存粗糙度
    outNormalRoughness = vec4(normalize(fragNormal), mat.roughnessFactor);
}