#version 460
layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outDenoisedSceneColor; // 现在输出的是最终合成画面！

layout(binding = 0) uniform sampler2D samplerNoisyGI;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;
// ★ 新增：接收直接光和贴图
layout(binding = 3) uniform sampler2D samplerDirectLight;
layout(binding = 4) uniform sampler2D samplerAlbedo;

layout(push_constant) uniform PushConstants {
    mat4 invViewProj;
    mat4 viewProj;
    vec3 cameraPos;
} pc;

const float sigmaSpatial = 4.0;
const float sigmaDepth = 0.1;
const float sigmaNormal = 32.0;


vec3 ReconstructWorldPos(vec2 uv, float depth) {
    vec4 ndc = vec4(uv * 2.0 - 1.0, depth, 1.0);
    vec4 worldPos = pc.invViewProj * ndc; 
    return worldPos.xyz / worldPos.w;
}

void main() {
    float exposure = 1.8; // 统一的曝光度
    float giMultiplier = 0.4; // GI 强度

    float centerDepth = texture(samplerDepth, inUV).r;
    
    // ★ 如果是天空，直接拿直接光（因为天空存在直接光里），加曝光并 ToneMapping
    if (centerDepth <= 0.0) {
        vec3 bgDirect = texture(samplerDirectLight, inUV).rgb;
        bgDirect = max(bgDirect, vec3(0.0)) * exposure;
        outDenoisedSceneColor = vec4(bgDirect, 1.0);
        return;
    }

    vec2 texSize = textureSize(samplerNoisyGI, 0);
    vec2 texelSize = 1.0 / texSize;

    vec3 centerWorldPos = ReconstructWorldPos(inUV, centerDepth);
    vec3 centerNormal = normalize(texture(samplerNormalRoughness, inUV).rgb + vec3(0.0001));
    vec3 centerGI = texture(samplerNoisyGI, inUV).rgb;

    vec3  sumGI = vec3(0.0);
    float sumWeight = 0.0;

    int radius = 2; // 降噪半径保持 2，配合 stepSize 足够了
    float stepSize = 3.0;

    for (int x = -radius; x <= radius; ++x) {
        for (int y = -radius; y <= radius; ++y) {
            vec2 offset = vec2(x, y) * stepSize;
            vec2 sampleUV = inUV + offset * texelSize;

            if (sampleUV.x < 0.0 || sampleUV.x > 1.0 || sampleUV.y < 0.0 || sampleUV.y > 1.0) continue;

            float sampleDepth = texture(samplerDepth, sampleUV).r;
            if (sampleDepth <= 0.0) continue;

            vec3 sampleGI = texture(samplerNoisyGI, sampleUV).rgb;
            vec3 sampleNormal = normalize(texture(samplerNormalRoughness, sampleUV).rgb + vec3(0.0001));
            vec3 sampleWorldPos = ReconstructWorldPos(sampleUV, sampleDepth);

            float distSq = dot(offset, offset);
            float wSpatial = exp(-distSq / (2.0 * sigmaSpatial * sigmaSpatial));

            vec3 dir = sampleWorldPos - centerWorldPos;
            float planeDist = abs(dot(centerNormal, dir));
            float wDepth = exp(-planeDist / sigmaDepth);

            float NdotN = max(dot(centerNormal, sampleNormal), 0.0);
            float wNormal = pow(NdotN, sigmaNormal);

            float weight = wSpatial * wDepth * wNormal;
            sumGI += sampleGI * weight;
            sumWeight += weight;
        }
    }

    // 1. 算出平滑后的纯净 GI
    vec3 finalGI = sumGI / max(sumWeight, 0.0001);

    // 2. 读取当前像素的材质和光照
    vec3 direct = texture(samplerDirectLight, inUV).rgb;
    vec3 albedo = texture(samplerAlbedo, inUV).rgb;

    // 3. 终极物理合成 + 曝光 + ACES
    vec3 sceneColor = direct + albedo * finalGI * giMultiplier;
    sceneColor = max(sceneColor, vec3(0.0)) * exposure;

    // 直接输出可以直接显示的完美图像！
    outDenoisedSceneColor = vec4(sceneColor, 1.0);
}