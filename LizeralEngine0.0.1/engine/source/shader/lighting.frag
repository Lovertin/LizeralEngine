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

layout(push_constant) uniform PushConstants {
    mat4 invViewProj;
    mat4 viewProj;
    vec3 cameraPos;
    float padding;
    uint frameIndex;
    uint padding2;
    uint64_t instanceDescAddr;
    vec3 lightDir;       
    float lightPadding;
    vec3 lightColor;
    float lightIntensity;
} pc;

struct RTInstanceDesc {
    uint64_t vertexBuffer;
    uint64_t indexBuffer;
    uint64_t materialBuffer;
    uint textureOffset;
    uint pad0, pad1, pad2;
};

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
    vec4 worldPos = pc.invViewProj * ndc; 
    return worldPos.xyz / worldPos.w;
}

layout(buffer_reference, scalar) readonly buffer VertexBuffer { float vData[]; };
layout(buffer_reference, scalar) readonly buffer IndexBuffer { uint iData[]; };
layout(buffer_reference, scalar) readonly buffer InstanceDescBuffer { RTInstanceDesc instances[]; };

vec3 TraceGlobalIlluminationRay(vec3 origin, vec3 direction) {
    rayQueryEXT rayQuery;
    uint rayFlags = gl_RayFlagsOpaqueEXT; // 不加 TerminateOnFirstHit，需要准确碰撞点
    rayQueryInitializeEXT(rayQuery, topLevelAS, rayFlags, 0xFF, origin, 0.01, direction, 1000.0);
    while(rayQueryProceedEXT(rayQuery)) {}

    if (rayQueryGetIntersectionTypeEXT(rayQuery, true) == gl_RayQueryCommittedIntersectionTriangleEXT) {
        // 1. 获取撞击信息
        uint instanceID = rayQueryGetIntersectionInstanceCustomIndexEXT(rayQuery, true);
        uint primID = rayQueryGetIntersectionPrimitiveIndexEXT(rayQuery, true);
        vec2 bary = rayQueryGetIntersectionBarycentricsEXT(rayQuery, true);

        // 2. 查台账
        InstanceDescBuffer instances = InstanceDescBuffer(pc.instanceDescAddr);
        RTInstanceDesc desc = instances.instances[instanceID];

        // 3. 查索引
        IndexBuffer indices = IndexBuffer(desc.indexBuffer);
        uint i0 = indices.iData[primID * 3 + 0];
        uint i1 = indices.iData[primID * 3 + 1];
        uint i2 = indices.iData[primID * 3 + 2];

        // 4. 查顶点计算 UV
        VertexBuffer verts = VertexBuffer(desc.vertexBuffer);
        vec2 uv0 = vec2(verts.vData[i0 * 8 + 6], verts.vData[i0 * 8 + 7]);
        vec2 uv1 = vec2(verts.vData[i1 * 8 + 6], verts.vData[i1 * 8 + 7]);
        vec2 uv2 = vec2(verts.vData[i2 * 8 + 6], verts.vData[i2 * 8 + 7]);

        vec2 hitUV = uv0 * (1.0 - bary.x - bary.y) + uv1 * bary.x + uv2 * bary.y;

        // 5. 采样无绑定贴图
        vec3 hitColor = texture(bindlessTextures[nonuniformEXT(desc.textureOffset)], hitUV).rgb;
        
        // ★ 核心除噪：去饱和度，压制采样到法线贴图带来的五彩斑斓的光污染
        float lum = dot(hitColor, vec3(0.299, 0.587, 0.114));
        hitColor = mix(hitColor, vec3(lum), 0.9);
        
        return hitColor; 
    }
    
    // 如果射向了虚空，返回极亮的天空颜色作为光照贡献
    return vec3(0.5, 0.7, 0.9) * 3.0; 
}

void main() {
    float depth = texture(samplerDepth, inUV).r;
    
    // =======================================================
    // 渲染背景天空
    // =======================================================
    if (depth <= 0.0) {
        vec3 bg = mix(vec3(0.1, 0.2, 0.3), vec3(0.01, 0.02, 0.05), inUV.y);
        outDirectLight = vec4(bg, 1.0);        // 背景色直接塞进光照通道
        outNoisyGI = vec4(0.0, 0.0, 0.0, 1.0); // 天空没有间接反弹光
        return;
    }

    vec4 albedoMetallic = texture(samplerAlbedoMetallic, inUV);
    vec4 normalRoughness = texture(samplerNormalRoughness, inUV);

    vec3 albedo = albedoMetallic.rgb; 
    float metallic = clamp(albedoMetallic.a, 0.0, 0.9);
    vec3 normal = normalize(normalRoughness.rgb + vec3(0.0001));
    float roughness = clamp(normalRoughness.a, 0.05, 1.0); 

    vec3 worldPos = ReconstructWorldPos(inUV, depth);
    vec3 viewDir = normalize(pc.cameraPos - worldPos);

    vec3 lightDir = normalize(pc.lightDir);
    vec3 lightColor = pc.lightColor * pc.lightIntensity;

    // =======================================================
    // 物理软阴影 (Ray Traced PCSS)
    // =======================================================
    uint baseSeed = uint(gl_FragCoord.x) * 1973u + uint(gl_FragCoord.y) * 9277u + pc.frameIndex * 26699u;
    float lightRadius = 0.05; 

    vec3 lightUp = abs(lightDir.y) > 0.99 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
    vec3 lightRight = normalize(cross(lightUp, lightDir));
    lightUp = cross(lightDir, lightRight);

    vec3 rayOrigin = worldPos + normal * 0.05; 
    float tMin = 0.001;     
    float tMax = 1000.0;    
    uint rayFlags = gl_RayFlagsTerminateOnFirstHitEXT | gl_RayFlagsOpaqueEXT | gl_RayFlagsSkipClosestHitShaderEXT;

    int NUM_SAMPLES = 8; 
    float shadowSum = 0.0; 

    float randomAngle = rand(baseSeed) * 2.0 * 3.14159265;
    const float GOLDEN_ANGLE = 2.39996323;

    for(int i = 0; i < NUM_SAMPLES; i++) {
        float r = sqrt((float(i) + 0.5) / float(NUM_SAMPLES)) * lightRadius; 
        float theta = float(i) * GOLDEN_ANGLE + randomAngle;
        vec2 diskOffset = vec2(r * cos(theta), r * sin(theta));
        vec3 jitteredRayDir = normalize(lightDir + lightRight * diskOffset.x + lightUp * diskOffset.y);

        rayQueryEXT rayQuery;
        rayQueryInitializeEXT(rayQuery, topLevelAS, rayFlags, 0xFF, rayOrigin, tMin, jitteredRayDir, tMax);
        while(rayQueryProceedEXT(rayQuery)) {}

        if (rayQueryGetIntersectionTypeEXT(rayQuery, true) == gl_RayQueryCommittedIntersectionNoneEXT) {
            shadowSum += 1.0; 
        }
    }

    float shadow = shadowSum / float(NUM_SAMPLES);

    // =======================================================
    // 拆分光照输出 1：带贴图颜色的纯净直接光
    // =======================================================
    float NdotL = max(dot(normal, lightDir), 0.0);
    vec3 halfVec = normalize(lightDir + viewDir);
    float NdotH = max(dot(normal, halfVec), 0.0);
    float shininess = exp2(10.0 * (1.0 - roughness) + 1.0);
    float specTerm = pow(NdotH, shininess) * (1.0 - roughness);
    
    vec3 specularColor = mix(vec3(1.0), albedo, metallic);
    
    vec3 pureDirectDiffuse = NdotL * lightColor; 
    vec3 pureDirectSpecular = specularColor * specTerm * lightColor; 
    
    // ★ 修复处：将 albedo 乘给了纯粹的光线！
    vec3 directLight = (albedo * pureDirectDiffuse * (1.0 - metallic) + pureDirectSpecular) * shadow;
    outDirectLight = vec4(directLight, 1.0);

    // =======================================================
    // 拆分光照输出 2：多采样结构化 GI (无贴图颜色)
    // =======================================================
    uint giSeed = baseSeed + 8888u; 
    vec3 w_tbn = normal;
    vec3 u_tbn = normalize(cross((abs(w_tbn.x) > 0.1 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0)), w_tbn));
    vec3 v_tbn = cross(w_tbn, u_tbn);
    mat3 tbn_hemisphere = mat3(u_tbn, v_tbn, w_tbn);

    int numGisamples = 4;
    vec3 giSum = vec3(0.0);
    const float phi = (1.0 + sqrt(5.0)) / 2.0; 

    for(int i=0; i<numGisamples; i++) {
        float u_vogel = mod(float(i) / phi, 1.0);
        float v_vogel = float(i) / float(numGisamples);
        
        float vogel_r = sqrt(u_vogel);
        float vogel_theta = v_vogel * 2.0 * 3.14159265 + randomAngle; 
        vec2 vogel_point = vec2(vogel_r * cos(vogel_theta), vogel_r * sin(vogel_theta));

        vec3 hemisphere_dir_tbn = vec3(vogel_point.x, vogel_point.y, sqrt(max(0.0, 1.0 - dot(vogel_point, vogel_point))));
        vec3 giDirection = normalize(tbn_hemisphere * hemisphere_dir_tbn);

        vec3 col = TraceGlobalIlluminationRay(worldPos + normal * 0.05, giDirection);
        
        // 防止异常极亮射线的萤火虫现象
        giSum += min(col, vec3(2.5)); 
    }
    
    vec3 bouncedColor = giSum / float(numGisamples);
    float giMultiplier = pc.lightIntensity * 0.25; 
    
    // 只输出物理反弹光度，不含任何木纹纹理
    vec3 indirectLight = bouncedColor * giMultiplier; 
    outNoisyGI = vec4(indirectLight, 1.0);
}