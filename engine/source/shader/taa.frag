#version 460
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require

layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor; 

layout(binding = 0) uniform sampler2D samplerDenoisedGI;
layout(binding = 1) uniform sampler2D samplerHistory;
layout(binding = 2) uniform sampler2D samplerVelocity;
layout(binding = 3) uniform sampler2D samplerAlbedo;    
layout(binding = 4) uniform sampler2D samplerDirectLight;

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

layout(push_constant) uniform PushConstants {
    FrameDataBuffer frameDataAddr;
    uint64_t instanceDescAddr;
    uint64_t pointLightsAddr;
    uint stepSize;
    uint pointLightCount;
} pc;

vec3 RGBToYCoCg(vec3 rgb) {
    return vec3(
         0.25 * rgb.r + 0.5 * rgb.g + 0.25 * rgb.b,
         0.5  * rgb.r               - 0.5  * rgb.b,
        -0.25 * rgb.r + 0.5 * rgb.g - 0.25 * rgb.b
    );
}

vec3 YCoCgToRGB(vec3 ycocg) {
    return vec3(
        ycocg.x + ycocg.y - ycocg.z,
        ycocg.x + ycocg.z,
        ycocg.x - ycocg.y - ycocg.z
    );
}

void main() {
    vec3 gi = texture(samplerDenoisedGI, inUV).rgb;
    vec3 albedo = texture(samplerAlbedo, inUV).rgb;
    vec3 direct = texture(samplerDirectLight, inUV).rgb;

    vec3 currentSceneColor = direct + albedo * gi;
    if (pc.frameDataAddr.frame.frameIndex < 3u) {
        outColor = vec4(currentSceneColor, 1.0);
        return;
    }
    
    vec2 velocity = texture(samplerVelocity, inUV).xy;
    velocity = clamp(velocity, vec2(-1.0), vec2(1.0));
    if (any(isnan(velocity)) || any(isinf(velocity))) {
        outColor = vec4(currentSceneColor, 1.0);
        return;
    }
    
    vec2 prevUV = inUV - velocity;
    float velocityLength = length(velocity);
    if (prevUV.x < 0.0 || prevUV.x > 1.0 || prevUV.y < 0.0 || prevUV.y > 1.0 || velocityLength > 0.2) {
        outColor = vec4(currentSceneColor, 1.0); 
        return;
    }

    vec4 historyColor = texture(samplerHistory, prevUV);
    if (any(isnan(historyColor.rgb)) || any(isinf(historyColor.rgb))) {
        outColor = vec4(currentSceneColor, 1.0);
        return;
    }

    vec2 texelSize = 1.0 / textureSize(samplerDenoisedGI, 0);
    vec3 m1 = vec3(0.0);
    vec3 m2 = vec3(0.0);

    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            vec3 neighborGI = texture(samplerDenoisedGI, inUV + vec2(x, y) * texelSize).rgb;
            vec3 neighborAlbedo = texture(samplerAlbedo, inUV + vec2(x, y) * texelSize).rgb;
            vec3 neighborDirect = texture(samplerDirectLight, inUV + vec2(x, y) * texelSize).rgb;
            
            vec3 nColor = neighborDirect + neighborAlbedo * neighborGI;
            vec3 ycocg = RGBToYCoCg(nColor);
            m1 += ycocg;
            m2 += ycocg * ycocg;
        }
    }

    m1 /= 9.0;
    m2 /= 9.0;

    vec3 mu = m1;
    vec3 sigma = sqrt(max(m2 - m1 * m1, vec3(0.0)));

    float gamma = 1.2; 
    vec3 colorMin = mu - gamma * sigma;
    vec3 colorMax = mu + gamma * sigma;

    vec3 historyYCoCg = RGBToYCoCg(historyColor.rgb);
    historyYCoCg = clamp(historyYCoCg, colorMin, colorMax);
    historyColor.rgb = YCoCgToRGB(historyYCoCg);

    float blendWeight = 0.12; 
    if (velocityLength > 0.001) {
        blendWeight = mix(0.12, 0.22, clamp(velocityLength * 40.0, 0.0, 1.0));
    }

    outColor = mix(historyColor, vec4(currentSceneColor, 1.0), blendWeight);
}
