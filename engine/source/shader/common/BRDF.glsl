// shader/Common/BRDF.glsl

#ifndef BRDF_GLSL
#define BRDF_GLSL

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
    float NDF = alpha2 / (PI * denom * denom); // 依赖 Math.glsl 中的 PI
    
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
    vec3 diffuse = kD * albedo / PI;

    return diffuse + specular;
}

#endif