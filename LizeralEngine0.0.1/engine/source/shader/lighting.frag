#version 460
// =======================================================
// ★ 1. 开启 Vulkan 光线查询 (Ray Query) 扩展
// =======================================================
#extension GL_EXT_ray_query : require

layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D samplerAlbedoMetallic;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;

// =======================================================
// ★ 2. 接收来自 C++ 端的场景大沙盘 (TLAS)
// 注意 binding = 4，与 C++ 中 BindAccelerationStructure 的绑定点一致
// =======================================================
layout(binding = 4) uniform accelerationStructureEXT topLevelAS;

layout(push_constant) uniform PushConstants {
    mat4 invViewProj;
    mat4 viewProj;
    vec3 cameraPos;
} pc;

// 坐标重构
vec3 ReconstructWorldPos(vec2 uv, float depth) {
    vec4 ndc = vec4(uv * 2.0 - 1.0, depth, 1.0);
    vec4 worldPos = pc.invViewProj * ndc; 
    return worldPos.xyz / worldPos.w;
}

// ACES 色调映射
vec3 ACESFilm(vec3 x) {
    float a = 2.51f; float b = 0.03f; float c = 2.43f; float d = 0.59f; float e = 0.14f;
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0, 1.0);
}

void main() {
    float depth = texture(samplerDepth, inUV).r;
    
    // 背景剔除
    if (depth <= 0.0) {
        vec3 bg = mix(vec3(0.1, 0.2, 0.3), vec3(0.01, 0.02, 0.05), inUV.y);
        outColor = vec4(bg, 1.0);
        return;
    }

    // 解码 G-Buffer
    vec4 albedoMetallic = texture(samplerAlbedoMetallic, inUV);
    vec4 normalRoughness = texture(samplerNormalRoughness, inUV);

    vec3 albedo = albedoMetallic.rgb; 
    float metallic = clamp(albedoMetallic.a, 0.0, 0.9);
    vec3 normal = normalize(normalRoughness.rgb + vec3(0.0001));
    float roughness = clamp(normalRoughness.a, 0.05, 1.0); 

    vec3 worldPos = ReconstructWorldPos(inUV, depth);
    vec3 viewDir = normalize(pc.cameraPos - worldPos);

    // 假设光源方向 (指向光源)
    vec3 lightDir = normalize(vec3(1.0, 2.0, 1.0));
    vec3 lightColor = vec3(2.0, 1.9, 1.8);

    // =======================================================
    // ★ 3. 硬件光线追踪：阴影计算
    // =======================================================
    float shadow = 1.0; // 默认 1.0 表示“被照亮”
    
    // 防阴影痤疮 (Shadow Acne)：把发射点沿着法线往外推一点点，防止光线打到自己身上
    vec3 rayOrigin = worldPos + normal * 0.05; 
    vec3 rayDir = lightDir; // 射线方向直接指向太阳
    float tMin = 0.001;     // 起点盲区
    float tMax = 1000.0;    // 射线最大射程 (太阳在无限远处)

    // 声明一个光线查询对象
    rayQueryEXT rayQuery;

    // 配置光线参数：
    // gl_RayFlagsTerminateOnFirstHitEXT : 极致优化！只要光线撞到了任何东西，立刻停止计算。
    // 因为对于阴影来说，被一个东西挡住和被一万个东西挡住，结果都是黑的！
    uint rayFlags = gl_RayFlagsTerminateOnFirstHitEXT | gl_RayFlagsOpaqueEXT | gl_RayFlagsSkipClosestHitShaderEXT;

    // 发射射线！
    rayQueryInitializeEXT(rayQuery, topLevelAS, rayFlags, 0xFF, rayOrigin, tMin, rayDir, tMax);

    // 让硬件在 BVH 树中步进求交
    while(rayQueryProceedEXT(rayQuery)) {}

    // 检查最后这根射线是不是“撞死”了
    if (rayQueryGetIntersectionTypeEXT(rayQuery, true) != gl_RayQueryCommittedIntersectionNoneEXT) {
        shadow = 0.0; // 撞到了遮挡物，这里是死黑的阴影！
    }

    // =======================================================
    // 4. Blinn-Phong 光照混合
    // =======================================================
    float NdotL = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = albedo * NdotL * lightColor;

    vec3 halfVec = normalize(lightDir + viewDir);
    float NdotH = max(dot(normal, halfVec), 0.0);
    float shininess = exp2(10.0 * (1.0 - roughness) + 1.0);
    float specTerm = pow(NdotH, shininess) * (1.0 - roughness);
    vec3 specularColor = mix(vec3(1.0), albedo, metallic);
    vec3 specular = specularColor * specTerm * lightColor;

    vec3 ambient = albedo * 0.15; // 环境光 (GI的下位替代)，不受直射光阴影影响

    // ★ 把光追算出来的 shadow 乘到漫反射和高光上！
    vec3 finalColor = ambient + (diffuse * (1.0 - metallic) + specular) * shadow;

    finalColor = ACESFilm(finalColor);
    outColor = vec4(finalColor, 1.0); 
}