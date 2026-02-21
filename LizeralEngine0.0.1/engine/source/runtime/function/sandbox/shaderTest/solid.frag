#version 330 core
out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoords;

// --- 摄像机与光照 ---
uniform vec3 u_camPos;
uniform vec3 u_lightDir;   // 指向光源的方向
uniform vec3 u_lightColor; // 光源颜色和强度

// --- PBR 材质参数 ---
uniform vec3  u_Albedo;
uniform float u_Metallic;
uniform float u_Roughness;
uniform float u_AO;

// 贴图采样器和开关宏
uniform sampler2D u_AlbedoMap;
uniform bool u_UseAlbedoMap;

uniform samplerCube u_IrradianceMap; // 漫反射环境光
uniform samplerCube u_PrefilterMap;  // 高光预滤波环境光
uniform bool u_UseIBL;               // IBL 开关

const float PI = 3.14159265359;

// 1. D: 法线分布函数 (GGX)
float DistributionGGX(vec3 N, vec3 H, float roughness) {
    float a = roughness * roughness; // 采用迪士尼原则，粗糙度平方视觉上更线性
    float a2 = a * a;
    float NdotH = max(dot(N, H), 0.0);
    float NdotH2 = NdotH * NdotH;

    float num = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return num / denom;
}

// 2. G: 几何函数 (Schlick-GGX)
float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = (roughness + 1.0);
    float k = (r * r) / 8.0;

    float num = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return num / denom;
}
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness) {
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);
    return ggx1 * ggx2;
}

// 3. F: 菲涅尔方程 (Fresnel-Schlick)
vec3 fresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

vec3 fresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness) {
    return F0 + (max(vec3(1.0 - roughness), F0) - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

void main() {

    // 动态决定使用纯色还是贴图
    vec3 albedoColor = u_Albedo;
    if (u_UseAlbedoMap) {
        // 从贴图中采样颜色。注意：PBR 要求 Albedo 必须在真实的线性空间中运算
        // 所以我们要用 pow 2.2 把 sRGB 的图片转换回线性空间
        albedoColor = pow(texture(u_AlbedoMap, TexCoords).rgb, vec3(2.2));
    }

    vec3 N = normalize(Normal);
    vec3 V = normalize(u_camPos - FragPos);

    // 基础反射率：非金属固定为 0.04，金属则使用反照率颜色
    vec3 F0 = vec3(0.04); 
    F0 = mix(F0, u_Albedo, u_Metallic);

    // 计算光照辐射率 (目前假设是平行光)
    vec3 L = normalize(u_lightDir);
    vec3 H = normalize(V + L); // 半程向量
    vec3 radiance = u_lightColor;

    // --- Cook-Torrance BRDF ---
    float NDF = DistributionGGX(N, H, u_Roughness);        
    float G   = GeometrySmith(N, V, L, u_Roughness);      
    vec3 F    = fresnelSchlick(max(dot(H, V), 0.0), F0);       

    vec3 numerator    = NDF * G * F;
    float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.0001; // 防止除以0
    vec3 specular     = numerator / denominator;

    // 能量守恒：反射的光 (kS) 和折射的光 (kD) 之和不能超过 1
    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS;
    kD *= 1.0 - u_Metallic; // 金属没有漫反射，全被吸收或反射了

    // 最终的输出光照
    float NdotL = max(dot(N, L), 0.0);        
    vec3 Lo = (kD * albedoColor / PI + specular) * radiance * NdotL;

    // 加上一点点环境光，防止完全死黑
    vec3 ambient;
    
    if (u_UseIBL) {
        // 1. 漫反射 IBL：根据法线 N 采样 Irradiance Map
        vec3 irradiance = texture(u_IrradianceMap, N).rgb;
        vec3 diffuse    = irradiance * albedoColor;

        // 2. 高光 IBL：根据反射向量 R 采样 Prefilter Map
        vec3 R = reflect(-V, N); // 视线打在表面后的反射方向
        // 假设预滤波贴图有 5 个 Mipmap 层级 (0.0 到 4.0)，根据粗糙度去查对应的模糊层级
        const float MAX_REFLECTION_LOD = 4.0;
        vec3 prefilteredColor = textureLod(u_PrefilterMap, R, u_Roughness * MAX_REFLECTION_LOD).rgb;

        // 环境光菲涅尔
        vec3 F_ambient = fresnelSchlickRoughness(max(dot(N, V), 0.0), F0, u_Roughness);
        
        vec3 kS_ambient = F_ambient;
        vec3 kD_ambient = 1.0 - kS_ambient;
        kD_ambient *= 1.0 - u_Metallic;

        // (注意：这里用 F_ambient 粗略近似了 BRDF 积分，为了快速出效果)
        vec3 specular_ambient = prefilteredColor * F_ambient;

        // 最终的 IBL 环境光
        ambient = (kD_ambient * diffuse + specular_ambient) * u_AO;
    } else {
        // 降级方案：如果没有 IBL 贴图，就给一点点固定的死光
        ambient = vec3(0.03) * albedoColor * u_AO;
    }
    
    vec3 color = ambient + Lo;

    // HDR 色调映射 (Tone Mapping) 和 Gamma 校正
    color = color / (color + vec3(1.0));
    color = pow(color, vec3(1.0/2.2)); 

    FragColor = vec4(color, 1.0);
}