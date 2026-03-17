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
    return vec3(0.5, 0.7, 0.9) * 0.2; 
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

    int NUM_SAMPLES = 4; 
    float shadowSum = 0.0; 

    float randomAngle = rand(baseSeed) * 2.0 * 3.14159265;
    const float GOLDEN_ANGLE = 2.39996323;
    float baseAngle = rand(baseSeed) * 2.0 * 3.14159265;

    for(int i = 0; i < NUM_SAMPLES; i++) {
        // Vogel Disk 算法：按斐波那契螺旋均匀且随机地散布在圆盘上
        float r = sqrt((float(i) + 0.5) / float(NUM_SAMPLES)) * lightRadius;
        float theta = float(i) * GOLDEN_ANGLE + baseAngle;
        
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
    shadow = pow(shadow, 2.0);

    float NdotL = max(dot(normal, lightDir), 0.0);
    float NdotV = max(dot(normal, viewDir), 0.0);
    vec3 H = normalize(lightDir + viewDir);
    float NdotH = max(dot(normal, H), 0.0);
    float VdotH = max(dot(viewDir, H), 0.0);

    // 1. 基础反射率 F0 (非金属一般为 0.04，金属则为自身颜色)
    vec3 F0 = mix(vec3(0.04), albedo, metallic);

    // 2. Fresnel 菲涅尔项 (Schlick近似) - 决定不同角度下的反光强度
    vec3 F = F0 + (1.0 - F0) * pow(clamp(1.0 - VdotH, 0.0, 1.0), 5.0);

    // 3. Normal Distribution 法线分布 (GGX) - 决定高光的集中度
    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
    float denom = (NdotH * NdotH * (alpha2 - 1.0) + 1.0);
    float NDF = alpha2 / (3.14159265 * denom * denom);

    // 4. Geometry 几何遮蔽 (Smith) - 决定微小凹凸的自阴影
    float k = (roughness + 1.0) * (roughness + 1.0) / 8.0;
    float gl = NdotL / (NdotL * (1.0 - k) + k);
    float gv = NdotV / (NdotV * (1.0 - k) + k);
    float G = gl * gv;

    // 5. 组合 BRDF
    vec3 numerator = NDF * G * F;
    float denominator = 4.0 * max(NdotV, 0.0) * max(NdotL, 0.0) + 0.0001;
    vec3 specular = numerator / denominator;

    // 漫反射 (能量守恒：非金属且没被反射走的光，才会被吸收并漫反射)
    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS;
    kD *= 1.0 - metallic;
    vec3 diffuse = kD * albedo / 3.14159265;

    // 最终直接光合成
    vec3 directLight = (diffuse + specular) * lightColor * NdotL * shadow;
    outDirectLight = vec4(directLight, 1.0);

    // =======================================================
    // 拆分光照输出 2：多采样结构化 GI (无贴图颜色)
    // =======================================================
    // 1. 构建切线空间矩阵 (TBN)
    vec3 w_tbn = normal;
    vec3 u_tbn = normalize(cross((abs(w_tbn.x) > 0.1 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0)), w_tbn));
    vec3 v_tbn = cross(w_tbn, u_tbn);
    mat3 tbn_hemisphere = mat3(u_tbn, v_tbn, w_tbn);

    int numGisamples = 4; // 采样数保持 4 不变
    vec3 giSum = vec3(0.0);

    // 2. 为了时空降噪，给每帧加上不同的偏移量 (Golden Ratio)
    // 这样同一个像素在不同帧会发射不同方向的射线，SVGF 就能把它们平均掉
    uint giSeed = baseSeed + pc.frameIndex * 161803u; 

    for(int i = 0; i < numGisamples; i++) {
        // 取两个 0~1 的随机数
        float r1 = rand(giSeed);
        float r2 = rand(giSeed);

        // ★ 数学魔法：将 2D 随机数转换为余弦加权的半球方向
        float r = sqrt(r1);
        float phi = 2.0 * 3.14159265 * r2;
        
        // 局部半球坐标
        vec3 localDir = vec3(
            r * cos(phi), 
            r * sin(phi), 
            sqrt(max(0.0, 1.0 - r1)) // Z 轴（法线方向）权重最大
        );

        // 转换到世界空间
        vec3 giDirection = normalize(tbn_hemisphere * localDir);

        // 发射光线
        vec3 col = TraceGlobalIlluminationRay(worldPos + normal * 0.05, giDirection);
        
        // 压制异常亮的“萤火虫”噪点
        giSum += min(col, vec3(5.0)); 
    }
    
    // 3. 计算最终颜色
    // 因为余弦采样的 PDF (概率密度) = cos(theta) / PI，它刚好和渲染方程的 cos 项抵消
    // 所以这里的公式变成了极简的：(Sum / N) * PI
    vec3 bouncedColor = (giSum / float(numGisamples)) * 3.14159265;
    
    // 你可以微调这个乘数来控制间接光的整体明暗
    float giMultiplier = pc.lightIntensity * 0.1; 
    
    vec3 indirectLight = bouncedColor * giMultiplier; 
    outNoisyGI = vec4(indirectLight, 1.0);
}