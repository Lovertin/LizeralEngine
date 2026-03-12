#version 460

layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D samplerAlbedoMetallic;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;

layout(push_constant) uniform PushConstants {
    mat4 invViewProj;
    mat4 viewProj;
    vec3 cameraPos;
} pc;

// =======================================================
// 坐标重构
// =======================================================
vec3 ReconstructWorldPos(vec2 uv, float depth) {
    vec4 ndc = vec4(uv * 2.0 - 1.0, depth, 1.0); 
    vec4 worldPos = pc.invViewProj * ndc; 
    return worldPos.xyz / worldPos.w; 
}

// ACES 色调映射
vec3 ACESFilm(vec3 x) {
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0, 1.0);
}

void main() {
    float depth = texture(samplerDepth, inUV).r;
    
    // 背景剔除 (如果是天空，直接输出天空蓝渐变)
    if (depth <= 0.0) {
        vec3 bg = mix(vec3(0.1, 0.2, 0.3), vec3(0.01, 0.02, 0.05), inUV.y);
        outColor = vec4(bg, 1.0);
        return;
    }

    vec4 albedoMetallic = texture(samplerAlbedoMetallic, inUV);
    vec4 normalRoughness = texture(samplerNormalRoughness, inUV);

    // 1. 数据解码 (线性化)
    vec3 albedo = albedoMetallic.rgb; 
    float metallic = clamp(albedoMetallic.a, 0.0, 0.9);
    
    // ★ 安全保护：强制归一化法线，防止模型导出错误导致法线为0引发 NaN 爆炸！
    vec3 normal = normalize(normalRoughness.rgb + vec3(0.0001)); 
    
    // 强制限制粗糙度，防止完全为 0 导致除以 0
    float roughness = clamp(normalRoughness.a, 0.05, 1.0); 

    vec3 worldPos = ReconstructWorldPos(inUV, depth);
    vec3 viewDir = normalize(pc.cameraPos - worldPos);

    // =======================================================
    // 2. 纯粹的 Blinn-Phong 光照 (降级版)
    // =======================================================
    // 假设有一盏主光源 (定向光)
    vec3 lightDir = normalize(vec3(1.0, 2.0, 1.0));
    vec3 lightColor = vec3(2.0, 1.9, 1.8);

    // 漫反射
    float NdotL = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = albedo * NdotL * lightColor;

    // 高光
    vec3 halfVec = normalize(lightDir + viewDir);
    float NdotH = max(dot(normal, halfVec), 0.0);
    float shininess = exp2(10.0 * (1.0 - roughness) + 1.0); 
    float specTerm = pow(NdotH, shininess) * (1.0 - roughness);
    vec3 specularColor = mix(vec3(1.0), albedo, metallic);
    vec3 specular = specularColor * specTerm * lightColor;

    // 极其基础的环境光 (代替 GI)
    vec3 ambient = albedo * 0.15; 

    // 合成颜色 (无反射)
    vec3 finalColor = ambient + diffuse * (1.0 - metallic) + specular;

    // =======================================================
    // 3. 输出
    // =======================================================
    finalColor = ACESFilm(finalColor);
    
    // 如果你的屏幕偏暗，可以解开这行 Gamma 校正
    // finalColor = pow(finalColor, vec3(1.0 / 2.2));

    outColor = vec4(finalColor, 1.0); 
}