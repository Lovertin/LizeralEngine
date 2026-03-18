#version 460
#extension GL_EXT_nonuniform_qualifier : require 
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require

layout(location = 0) in vec3 fragNormal;
layout(location = 1) in vec2 fragUV;
layout(location = 2) flat in uint fragMatID; // ★ 材质 ID
layout(location = 3) in vec4 fragCurrentPosClip;
layout(location = 4) in vec4 fragPrevPosClip;
layout(location = 5) in vec3 fragWorldPos;   // ★ 世界坐标

layout(location = 2) out vec2 outVelocity;
layout(location = 0) out vec4 outAlbedoMetallic; // RT0
layout(location = 1) out vec4 outNormalRoughness; // RT1

layout(binding = 0) uniform sampler2D GlobalTextures[1024]; 

// ★ 全新升级的 PBR 材质结构
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

vec3 PerturbNormal(vec3 worldNormal, vec3 worldPos, vec2 uv, vec3 normalMapSample) {
    // 对世界坐标和 UV 进行屏幕空间求导 (算出像素间的变化率)
    vec3 q1  = dFdx(worldPos);
    vec3 q2  = dFdy(worldPos);
    vec2 st1 = dFdx(uv);
    vec2 st2 = dFdy(uv);

    vec3 N   = normalize(worldNormal);
    // 叉乘推导出完美的切线 T 和副切线 B
    vec3 T   = normalize(q1 * st2.t - q2 * st1.t);
    vec3 B   = -normalize(cross(N, T));
    mat3 TBN = mat3(T, B, N);

    // 将 0~1 的法线贴图解压为 -1~1，并转换到世界空间
    vec3 tNormal = normalMapSample * 2.0 - 1.0;
    return normalize(TBN * tNormal);
}

void main() {
    MaterialBuffer matBuf = MaterialBuffer(pc.matBuf);
    uint localMatID = fragMatID - pc.textureOffset; 
    Material mat = matBuf.m[localMatID];

    // ====================================================
    // 1. 解析 BaseColor (Albedo)
    // ====================================================
    vec4 albedo = mat.baseColorFactor;
    if (mat.albedoTex >= 0) {
        albedo *= texture(GlobalTextures[nonuniformEXT(mat.albedoTex)], fragUV);
    }

    // ====================================================
    // 2. 解析 ORM (Occlusion, Roughness, Metallic)
    // ====================================================
    float ao = 1.0;
    float roughness = mat.roughnessFactor;
    float metallic = mat.metallicFactor;
    
    if (mat.ormTex >= 0) {
        vec3 orm = texture(GlobalTextures[nonuniformEXT(mat.ormTex)], fragUV).rgb;
        ao = orm.r;              // R 通道存 AO
        roughness *= orm.g;      // G 通道存 Roughness
        metallic *= orm.b;       // B 通道存 Metallic
    }

    // ====================================================
    // 3. 解析 Normal Map
    // ====================================================
    vec3 finalNormal = normalize(fragNormal);
    if (mat.normalTex >= 0) {
        vec3 normalMapSample = texture(GlobalTextures[nonuniformEXT(mat.normalTex)], fragUV).rgb;
        finalNormal = PerturbNormal(finalNormal, fragWorldPos, fragUV, normalMapSample);
    }

    // ====================================================
    // 4. 计算运动向量 (TAA)
    // ====================================================
    vec2 currentNDC = fragCurrentPosClip.xy / fragCurrentPosClip.w;
    vec2 prevNDC = fragPrevPosClip.xy / fragPrevPosClip.w;
    vec2 currentUV = currentNDC * 0.5 + 0.5;
    vec2 prevUV = prevNDC * 0.5 + 0.5;
    outVelocity = currentUV - prevUV;

    // ====================================================
    // 5. 组装 G-Buffer！
    // ====================================================
    // 注意：我们将 AO 预乘进了 Albedo 中，这是一种针对间接光渲染极其高效的 Hack 做法
    // 在真实 3A 中，AO 通常单独存，然后在 Lighting Pass 仅乘给 GI
    vec3 finalAlbedo = albedo.rgb * ao;

    // RT0: RGB 存颜色，Alpha 存金属度
    outAlbedoMetallic = vec4(finalAlbedo, metallic);

    // RT1: RGB 存带细节的法线，Alpha 存粗糙度
    outNormalRoughness = vec4(finalNormal, roughness);
}