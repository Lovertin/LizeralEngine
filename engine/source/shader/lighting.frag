#version 460
#extension GL_EXT_ray_query : require
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_nonuniform_qualifier : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require

layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outDirectLight;
layout(location = 1) out vec4 outNoisyGI;

layout(binding = 0) uniform sampler2D samplerAlbedoMetallic;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;
layout(binding = 4) uniform accelerationStructureEXT topLevelAS;
layout(set = 1, binding = 0) uniform sampler2D bindlessTextures[];

// ==========================================================
// 1. 全局与点光源的 BDA 数据结构
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

layout(buffer_reference, scalar, buffer_reference_align = 16) readonly buffer FrameDataBuffer {
    GlobalFrameData frame;
};

// ★ 新增：点光源结构与 Buffer 指针
struct GPUPointLight {
    vec4 posAndRadius;      // xyz: 位置, w: 半径
    vec4 colorAndIntensity; // rgb: 颜色, a: 强度
};
layout(buffer_reference, scalar, buffer_reference_align = 16) readonly buffer PointLightBuffer { 
    GPUPointLight lights[]; 
};

// ==========================================================
// 2. 更新版的 Push Constants
// ==========================================================
layout(push_constant) uniform PushConstants {
    FrameDataBuffer frameData;
    uint64_t instanceDescAddr;       
    uint64_t pointLightsAddr;  // ★ 新增：点光源数组首地址
    uint stepSize;             
    uint pointLightCount;      // ★ 新增：当前帧激活的点光源数量
} pc;

struct RTInstanceDesc {
    uint64_t vertexBuffer;
    uint64_t indexBuffer;
    uint64_t materialBuffer;
    uint textureOffset;
    uint pad0, pad1, pad2;
};

layout(buffer_reference, scalar) readonly buffer VertexBuffer { float vData[]; };
layout(buffer_reference, scalar) readonly buffer IndexBuffer { uint iData[]; };
layout(buffer_reference, scalar) readonly buffer InstanceDescBuffer { RTInstanceDesc instances[]; };

// --- 随机数与工具函数 ---
uint pcg_hash(uint seed) {
    uint state = seed * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

float rand(inout uint seed) {
    seed = pcg_hash(seed);
    return float(seed & 0x00FFFFFFu) / 16777216.0;
}

vec3 ReconstructWorldPos(vec2 uv, float depth) {
    vec4 ndc = vec4(uv * 2.0 - 1.0, depth, 1.0);
    vec4 worldPos = pc.frameData.frame.invViewProj * ndc; 
    return worldPos.xyz / worldPos.w;
}

// --- 光追环境光遮蔽/反射 ---
vec3 TraceGlobalIlluminationRay(vec3 origin, vec3 direction) {
    rayQueryEXT rayQuery;
    uint rayFlags = gl_RayFlagsOpaqueEXT; 
    rayQueryInitializeEXT(rayQuery, topLevelAS, rayFlags, 0xFF, origin, 0.01, direction, 1000.0);
    while(rayQueryProceedEXT(rayQuery)) {}

    if (rayQueryGetIntersectionTypeEXT(rayQuery, true) == gl_RayQueryCommittedIntersectionTriangleEXT) {
        uint instanceID = rayQueryGetIntersectionInstanceCustomIndexEXT(rayQuery, true);
        uint primID = rayQueryGetIntersectionPrimitiveIndexEXT(rayQuery, true);
        vec2 bary = rayQueryGetIntersectionBarycentricsEXT(rayQuery, true);

        InstanceDescBuffer instances = InstanceDescBuffer(pc.instanceDescAddr);
        RTInstanceDesc desc = instances.instances[instanceID];
        IndexBuffer indices = IndexBuffer(desc.indexBuffer);
        
        uint i0 = indices.iData[primID * 3 + 0];
        uint i1 = indices.iData[primID * 3 + 1];
        uint i2 = indices.iData[primID * 3 + 2];

        VertexBuffer verts = VertexBuffer(desc.vertexBuffer);
        vec2 uv0 = vec2(verts.vData[i0 * 8 + 6], verts.vData[i0 * 8 + 7]);
        vec2 uv1 = vec2(verts.vData[i1 * 8 + 6], verts.vData[i1 * 8 + 7]);
        vec2 uv2 = vec2(verts.vData[i2 * 8 + 6], verts.vData[i2 * 8 + 7]);
        
        vec2 hitUV = uv0 * (1.0 - bary.x - bary.y) + uv1 * bary.x + uv2 * bary.y;
        vec3 hitColor = texture(bindlessTextures[nonuniformEXT(desc.textureOffset)], hitUV).rgb;
        
        float lum = dot(hitColor, vec3(0.299, 0.587, 0.114));
        hitColor = mix(hitColor, vec3(lum), 0.9);
        return hitColor;
    }
    return vec3(0.5, 0.7, 0.9) * 0.2;
}

// ==========================================================
// ★ 提取出来的标准 PBR (BRDF) 计算函数
// ==========================================================
vec3 CalculateBRDF(vec3 normal, vec3 viewDir, vec3 lightDir, vec3 albedo, float metallic, float roughness, out float NdotL) {
    vec3 H = normalize(lightDir + viewDir);
    NdotL = max(dot(normal, lightDir), 0.0);
    float NdotV = max(dot(normal, viewDir), 0.0);
    float NdotH = max(dot(normal, H), 0.0);
    float VdotH = max(dot(viewDir, H), 0.0);

    vec3 F0 = mix(vec3(0.04), albedo, metallic);
    vec3 F = F0 + (1.0 - F0) * pow(clamp(1.0 - VdotH, 0.0, 1.0), 5.0);

    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
    float denom = (NdotH * NdotH * (alpha2 - 1.0) + 1.0);
    float NDF = alpha2 / (3.14159265 * denom * denom);

    float k = (roughness + 1.0) * (roughness + 1.0) / 8.0;
    float gl = NdotL / (NdotL * (1.0 - k) + k);
    float gv = NdotV / (NdotV * (1.0 - k) + k);
    float G = gl * gv;

    vec3 numerator = NDF * G * F;
    float denominator = 4.0 * max(NdotV, 0.0) * max(NdotL, 0.0) + 0.0001;
    vec3 specular = numerator / denominator;

    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS;
    kD *= 1.0 - metallic;
    vec3 diffuse = kD * albedo / 3.14159265;

    return diffuse + specular;
}

void main() {
    float depth = texture(samplerDepth, inUV).r;
    if (depth <= 0.0) {
        vec3 bg = mix(vec3(0.1, 0.2, 0.3), vec3(0.01, 0.02, 0.05), inUV.y);
        outDirectLight = vec4(bg, 1.0);        
        outNoisyGI = vec4(0.0, 0.0, 0.0, 1.0); 
        return;
    }

    vec4 albedoMetallic = texture(samplerAlbedoMetallic, inUV);
    vec4 normalRoughness = texture(samplerNormalRoughness, inUV);

    vec3 albedo = albedoMetallic.rgb; 
    float metallic = clamp(albedoMetallic.a, 0.0, 0.9);
    vec3 normal = normalize(normalRoughness.rgb + vec3(0.0001));
    float roughness = clamp(normalRoughness.a, 0.05, 1.0); 

    vec3 worldPos = ReconstructWorldPos(inUV, depth);
    vec3 viewDir = normalize(pc.frameData.frame.cameraPos - worldPos);
    
    uint baseSeed = uint(gl_FragCoord.x) * 1973u + uint(gl_FragCoord.y) * 9277u + pc.frameData.frame.frameIndex * 26699u;
    vec3 rayOrigin = worldPos + normal * 0.05;
    float tMin = 0.001;     
    uint rayFlags = gl_RayFlagsTerminateOnFirstHitEXT | gl_RayFlagsOpaqueEXT | gl_RayFlagsSkipClosestHitShaderEXT;

    vec3 finalDirectLight = vec3(0.0);

    // ====================================================
    // 1. 处理平行光 (太阳光)
    // ====================================================
    vec3 mainLightDir = normalize(pc.frameData.frame.lightDir);
    float lightRadius = 0.05; 
    vec3 lightUp = abs(mainLightDir.y) > 0.99 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
    vec3 lightRight = normalize(cross(lightUp, mainLightDir));
    lightUp = cross(mainLightDir, lightRight);

    int NUM_SAMPLES = 4;
    float shadowSum = 0.0; 
    float baseAngle = rand(baseSeed) * 2.0 * 3.14159265;
    
    for(int i = 0; i < NUM_SAMPLES; i++) {
        float r = sqrt((float(i) + 0.5) / float(NUM_SAMPLES)) * lightRadius;
        float theta = float(i) * 2.39996323 + baseAngle;
        vec2 diskOffset = vec2(r * cos(theta), r * sin(theta));
        vec3 jitteredRayDir = normalize(mainLightDir + lightRight * diskOffset.x + lightUp * diskOffset.y);

        rayQueryEXT rayQuery;
        rayQueryInitializeEXT(rayQuery, topLevelAS, rayFlags, 0xFF, rayOrigin, tMin, jitteredRayDir, 1000.0);
        while(rayQueryProceedEXT(rayQuery)) {}
        if (rayQueryGetIntersectionTypeEXT(rayQuery, true) == gl_RayQueryCommittedIntersectionNoneEXT) {
            shadowSum += 1.0;
        }
    }
    float shadow = pow(shadowSum / float(NUM_SAMPLES), 2.0);

    float NdotL_main;
    vec3 mainBrdf = CalculateBRDF(normal, viewDir, mainLightDir, albedo, metallic, roughness, NdotL_main);
    finalDirectLight += mainBrdf * pc.frameData.frame.lightColor * pc.frameData.frame.lightIntensity * NdotL_main * shadow;

    // ====================================================
    // ★ 2. 处理动态点光源
    // ====================================================
    if (pc.pointLightsAddr != 0 && pc.pointLightCount > 0) {
        PointLightBuffer plBuf = PointLightBuffer(pc.pointLightsAddr);
        
        for (uint i = 0; i < pc.pointLightCount; i++) {
            GPUPointLight pl = plBuf.lights[i];
            vec3 lightPos = pl.posAndRadius.xyz;
            float radius = pl.posAndRadius.w;
            vec3 lightColor = pl.colorAndIntensity.rgb;
            float intensity = pl.colorAndIntensity.a;

            vec3 L = lightPos - worldPos;
            float distance = length(L);

            // 如果像素超出了点光源的包围球，直接跳过
            if (distance < radius) {
                L /= distance; // 归一化光照方向

                // UE4 物理平方反比衰减公式 (Inverse Square Falloff)
                float d_over_r = distance / radius;
                float d_over_r_4 = d_over_r * d_over_r * d_over_r * d_over_r;
                float falloff = clamp(1.0 - d_over_r_4, 0.0, 1.0);
                falloff = (falloff * falloff) / (distance * distance + 1.0);

                // 发射一根光线检测点光源遮挡 (硬阴影性能更好)
                float plShadow = 1.0;
                rayQueryEXT rq;
                rayQueryInitializeEXT(rq, topLevelAS, rayFlags, 0xFF, rayOrigin, tMin, L, distance);
                while(rayQueryProceedEXT(rq)) {}
                if (rayQueryGetIntersectionTypeEXT(rq, true) != gl_RayQueryCommittedIntersectionNoneEXT) {
                    plShadow = 0.0; // 被挡住了
                }

                float NdotL_pl;
                vec3 plBrdf = CalculateBRDF(normal, viewDir, L, albedo, metallic, roughness, NdotL_pl);
                
                // 累加点光源结果
                finalDirectLight += plBrdf * lightColor * intensity * NdotL_pl * falloff * plShadow;
            }
        }
    }

    outDirectLight = vec4(finalDirectLight, 1.0);

    // ====================================================
    // 3. GI 采样 (间接光不变)
    // ====================================================
    vec3 w_tbn = normal;
    vec3 u_tbn = normalize(cross((abs(w_tbn.x) > 0.1 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0)), w_tbn));
    vec3 v_tbn = cross(w_tbn, u_tbn);
    mat3 tbn_hemisphere = mat3(u_tbn, v_tbn, w_tbn);

    int numGisamples = 4; 
    vec3 giSum = vec3(0.0);
    uint giSeed = baseSeed + pc.frameData.frame.frameIndex * 161803u; 

    for(int i = 0; i < numGisamples; i++) {
        float r1 = rand(giSeed);
        float r2 = rand(giSeed);
        float r = sqrt(r1);
        float phi = 2.0 * 3.14159265 * r2;
        vec3 localDir = vec3(r * cos(phi), r * sin(phi), sqrt(max(0.0, 1.0 - r1)));
        vec3 giDirection = normalize(tbn_hemisphere * localDir);
        vec3 col = TraceGlobalIlluminationRay(worldPos + normal * 0.05, giDirection);
        giSum += min(col, vec3(5.0));
    }
    
    vec3 bouncedColor = (giSum / float(numGisamples)) * 3.14159265;
    float giMultiplier = pc.frameData.frame.lightIntensity * 0.04; 
    vec3 indirectLight = bouncedColor * giMultiplier;
    outNoisyGI = vec4(indirectLight, 1.0);
}