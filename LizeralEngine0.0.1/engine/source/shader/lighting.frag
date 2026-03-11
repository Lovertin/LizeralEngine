#version 460

layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor;

// 接收我们刚刚画好的 G-Buffer
layout(binding = 0) uniform sampler2D samplerAlbedoMetallic;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;

void main() {
    // 1. 从 G-Buffer 采样数据
    vec4 albedoMetallic = texture(samplerAlbedoMetallic, inUV);
    vec4 normalRoughness = texture(samplerNormalRoughness, inUV);
    float depth = texture(samplerDepth, inUV).r;

    // 2. 解码
    vec3 albedo = albedoMetallic.rgb;
    float metallic = albedoMetallic.a;
    vec3 normal = normalRoughness.rgb;
    float roughness = normalRoughness.a;

    // 如果是背景（深度为 1.0），直接输出背景色
    if (depth >= 1.0) {
        outColor = vec4(0.2, 0.2, 0.2, 1.0); // 深灰色背景
        return;
    }

    // 3. 极其纯粹的光照计算！(未来你可以把 Lumen 算法写在这里)
    vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
    float diff = max(dot(normal, lightDir), 0.0);
    float ambient = mix(0.4, 0.1, roughness); 

    outColor = vec4(albedo * (diff + ambient), 1.0);
}