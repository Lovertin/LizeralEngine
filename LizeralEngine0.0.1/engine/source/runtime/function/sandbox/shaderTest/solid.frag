#version 330 core
out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoords;
in vec4 FragPosViewSpace; // 接收观察空间坐标

// --- 摄像机与光照 ---
uniform vec3 u_camPos;
uniform vec3 u_lightDir;
uniform vec3 u_lightColor;

// --- PBR 材质参数 ---
uniform vec3  u_Albedo;
uniform float u_Metallic;
uniform float u_Roughness;
uniform float u_AO;

uniform sampler2D u_AlbedoMap;
uniform bool u_UseAlbedoMap;

uniform samplerCube u_IrradianceMap;
uniform samplerCube u_PrefilterMap; 
uniform bool u_UseIBL;

uniform mat4 u_lightSpaceMatrices[4];     // 级联矩阵数组
uniform float u_cascadePlaneDistances[4]; // 级联分割深度
uniform int u_cascadeCount;               // 级联层数
uniform sampler2DArray shadowMap;         // 【注意】：变成了 2DArray 采样器！
uniform mat4 view;                        // 相机视图矩阵

const float PI = 3.14159265359;

const vec2 poissonDisk[16] = vec2[](
    vec2( -0.94201624, -0.39906216 ), vec2(  0.94558609, -0.76890725 ),
    vec2( -0.094184101, -0.92938870), vec2(  0.34495938,  0.29387760 ),
    vec2( -0.91588401,  0.45771432 ), vec2( -0.81544232, -0.87912464 ),
    vec2( -0.38277543,  0.27676845 ), vec2(  0.97484398,  0.75648379 ),
    vec2(  0.44323325, -0.97511554 ), vec2(  0.53742981, -0.47373420 ),
    vec2( -0.26496911, -0.41893023 ), vec2(  0.79197514,  0.19090188 ),
    vec2( -0.24188840,  0.99706507 ), vec2( -0.81409955,  0.91437590 ),
    vec2(  0.19984126,  0.78641367 ), vec2(  0.14383161, -0.14100790 )
);

// ... (DistributionGGX, GeometrySmith, fresnelSchlick 等函数保持不变) ...
float DistributionGGX(vec3 N, vec3 H, float roughness) {
    float a = roughness * roughness;
    float a2 = a * a;
    float NdotH = max(dot(N, H), 0.0);
    float NdotH2 = NdotH * NdotH;
    float num = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
    return num / denom;
}

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

vec3 fresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

vec3 fresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness) {
    return F0 + (max(vec3(1.0 - roughness), F0) - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

float ShadowCalculation(vec3 fragPosWorld, vec3 N, vec3 L) {
    // 1. 获取像素离摄像机的距离
    float depthValue = abs(FragPosViewSpace.z);

    // 2. 智能选层
    int layer = -1;
    for (int i = 0; i < u_cascadeCount; ++i) {
        if (depthValue < u_cascadePlaneDistances[i]) {
            layer = i;
            break;
        }
    }
    if (layer == -1) {
        layer = u_cascadeCount - 1;
    }

    float normalOffsetScale = 0.05; // 基础偏移量，可根据你的场景大小微调
    if (layer == 1) normalOffsetScale *= 1.5;
    else if (layer == 2) normalOffsetScale *= 2.0;
    else if (layer == 3) normalOffsetScale *= 3.0;

    // 光线越是倾斜照射表面，需要的法线偏移越大
    vec3 offsetPos = fragPosWorld + N * normalOffsetScale * (1.0 - dot(N, L));

    // 3. 用偏移后的坐标拿到该层对应的太阳矩阵
    vec4 fragPosLightSpace = u_lightSpaceMatrices[layer] * vec4(offsetPos, 1.0);
    
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
    projCoords = projCoords * 0.5 + 0.5;

    if(projCoords.z > 1.0)
        return 0.0;

    float currentDepth = projCoords.z;

    float bias = 0.0005;
    if (layer == 1) bias *= 1.2;
    else if (layer == 2) bias *= 1.5;
    else if (layer == 3) bias *= 2.0;

    // 4. 泊松圆盘 PCF
    float shadow = 0.0;
    vec2 texelSize = 1.0 / vec2(textureSize(shadowMap, 0).xy); 
    
    float filterRadius = 2.0; 
    if (layer == 1) filterRadius *= 0.5;  
    else if (layer == 2) filterRadius *= 0.25; 
    else if (layer == 3) filterRadius *= 0.1;

    for(int i = 0; i < 16; i++) {
        float pcfDepth = texture(shadowMap, vec3(projCoords.xy + poissonDisk[i] * texelSize * filterRadius, float(layer))).r; 
        shadow += currentDepth - bias > pcfDepth ? 1.0 : 0.0;        
    }
    shadow /= 16.0;
    
    return shadow;
}

void main() {
    vec3 albedoColor = u_Albedo;
    if (u_UseAlbedoMap) {
        albedoColor = pow(texture(u_AlbedoMap, TexCoords).rgb, vec3(2.2));
    }

    vec3 N = normalize(Normal);
    vec3 V = normalize(u_camPos - FragPos);

    vec3 F0 = vec3(0.04); 
    F0 = mix(F0, albedoColor, u_Metallic);

    vec3 L = normalize(u_lightDir);
    vec3 H = normalize(V + L); 
    vec3 radiance = u_lightColor;

    float NDF = DistributionGGX(N, H, u_Roughness);        
    float G   = GeometrySmith(N, V, L, u_Roughness);      
    vec3 F    = fresnelSchlick(max(dot(H, V), 0.0), F0);       

    vec3 numerator    = NDF * G * F;
    float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.0001; 
    vec3 specular     = numerator / denominator;

    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS;
    kD *= 1.0 - u_Metallic; 

    // 【修改】：只需传入世界坐标 FragPos
    float shadow = ShadowCalculation(FragPos, N, L);

    float NdotL = max(dot(N, L), 0.0);        
    vec3 Lo = (1.0 - shadow) * (kD * albedoColor / PI + specular) * radiance * NdotL;

    vec3 ambient;
    if (u_UseIBL) {
        vec3 irradiance = texture(u_IrradianceMap, N).rgb;
        vec3 diffuse    = irradiance * albedoColor;

        vec3 R = reflect(-V, N); 
        const float MAX_REFLECTION_LOD = 4.0;
        vec3 prefilteredColor = textureLod(u_PrefilterMap, R, u_Roughness * MAX_REFLECTION_LOD).rgb;

        vec3 F_ambient = fresnelSchlickRoughness(max(dot(N, V), 0.0), F0, u_Roughness);
        
        vec3 kS_ambient = F_ambient;
        vec3 kD_ambient = 1.0 - kS_ambient;
        kD_ambient *= 1.0 - u_Metallic;

        vec3 specular_ambient = prefilteredColor * F_ambient;
        ambient = (kD_ambient * diffuse + specular_ambient) * u_AO;
    } else {
        ambient = vec3(0.03) * albedoColor * u_AO;
    }
    
    vec3 color = ambient + Lo;
    color = color / (color + vec3(1.0));
    color = pow(color, vec3(1.0/2.2)); 

    FragColor = vec4(color, 1.0);
}