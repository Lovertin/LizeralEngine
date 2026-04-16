#version 460
#extension GL_EXT_nonuniform_qualifier : require
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require

layout(location = 0) in vec3 fragNormal;
layout(location = 1) in vec2 fragUV;
layout(location = 2) flat in uint fragMatID;
layout(location = 3) in vec4 fragCurrentPosClip;
layout(location = 4) in vec4 fragPrevPosClip;
layout(location = 5) in vec3 fragWorldPos;

layout(location = 0) out vec4 outColor;

layout(set = 0, binding = 0) uniform sampler2D samplerSceneDepth;
layout(set = 1, binding = 0) uniform sampler2D GlobalTextures[1024];

struct Material {
    vec4 baseColorFactor;
    float metallicFactor;
    float roughnessFactor;
    int albedoTex;
    int normalTex;
    int ormTex;
    int emissiveTex;
    int alphaMode;
    float alphaCutoff;
    int pad0;
    int pad1;
};

layout(buffer_reference, scalar, buffer_reference_align = 4) readonly buffer MaterialBuffer { Material m[]; };

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

layout(buffer_reference, scalar, buffer_reference_align = 16) readonly buffer InstanceDataBuffer {
    mat4 model;
    mat4 prevModel;
};

layout(push_constant) uniform PushConstants {
    FrameDataBuffer frameDataAddr;
    InstanceDataBuffer instanceDataAddr;
    uint64_t vBuf;
    uint64_t mBuf;
    uint64_t iBuf;
    uint64_t bBuf;
    uint64_t matBuf;
    uint totalMeshlets;
    uint textureOffset;
} pc;

vec3 PerturbNormal(vec3 worldNormal, vec3 worldPos, vec2 uv, vec3 normalMapSample) {
    vec3 q1 = dFdx(worldPos);
    vec3 q2 = dFdy(worldPos);
    vec2 st1 = dFdx(uv);
    vec2 st2 = dFdy(uv);

    vec3 N = normalize(worldNormal);
    vec3 T = normalize(q1 * st2.t - q2 * st1.t);
    vec3 B = -normalize(cross(N, T));
    mat3 TBN = mat3(T, B, N);

    vec3 tNormal = normalMapSample * 2.0 - 1.0;
    return normalize(TBN * tNormal);
}

vec3 FresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

void main() {
    vec2 screenUV = gl_FragCoord.xy / vec2(textureSize(samplerSceneDepth, 0));
    float sceneDepth = texture(samplerSceneDepth, screenUV).r;
    if (gl_FragCoord.z < sceneDepth - 1e-5) {
        discard;
    }

    MaterialBuffer matBuf = MaterialBuffer(pc.matBuf);
    Material mat = matBuf.m[fragMatID];

    if (mat.alphaMode != 2) {
        discard;
    }

    vec4 albedo = mat.baseColorFactor;
    if (mat.albedoTex >= 0) {
        albedo *= texture(GlobalTextures[nonuniformEXT(mat.albedoTex)], fragUV);
    }

    if (albedo.a <= 0.001) {
        discard;
    }

    float roughness = clamp(mat.roughnessFactor, 0.05, 1.0);
    float metallic = clamp(mat.metallicFactor, 0.0, 1.0);
    if (mat.ormTex >= 0) {
        vec3 orm = texture(GlobalTextures[nonuniformEXT(mat.ormTex)], fragUV).rgb;
        roughness *= orm.g;
        metallic *= orm.b;
    }
    roughness = clamp(roughness, 0.05, 1.0);

    vec3 N = normalize(fragNormal);
    if (mat.normalTex >= 0) {
        vec3 normalMapSample = texture(GlobalTextures[nonuniformEXT(mat.normalTex)], fragUV).rgb;
        N = PerturbNormal(N, fragWorldPos, fragUV, normalMapSample);
    }

    vec3 V = normalize(pc.frameDataAddr.frame.cameraPos - fragWorldPos);
    vec3 L = normalize(pc.frameDataAddr.frame.lightDir);
    vec3 H = normalize(V + L);
    float NdotL = max(dot(N, L), 0.0);
    float NdotV = max(dot(N, V), 0.0);
    float NdotH = max(dot(N, H), 0.0);
    float VdotH = max(dot(V, H), 0.0);

    vec3 baseColor = albedo.rgb;
    vec3 F0 = mix(vec3(0.04), baseColor, metallic);
    vec3 F = FresnelSchlick(VdotH, F0);

    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
    float denom = (NdotH * NdotH * (alpha2 - 1.0) + 1.0);
    float D = alpha2 / max(3.14159265 * denom * denom, 1e-4);

    float k = (roughness + 1.0);
    k = (k * k) / 8.0;
    float Gv = NdotV / max(NdotV * (1.0 - k) + k, 1e-4);
    float Gl = NdotL / max(NdotL * (1.0 - k) + k, 1e-4);
    float G = Gv * Gl;

    vec3 specular = (D * G * F) / max(4.0 * NdotV * NdotL, 1e-4);
    vec3 kS = F;
    vec3 kD = (vec3(1.0) - kS) * (1.0 - metallic);
    vec3 diffuse = kD * baseColor / 3.14159265;

    vec3 direct = (diffuse + specular) * pc.frameDataAddr.frame.lightColor * pc.frameDataAddr.frame.lightIntensity * NdotL;
    vec3 ambient = baseColor * 0.03;
    vec3 color = ambient + direct;

    if (mat.emissiveTex >= 0) {
        color += texture(GlobalTextures[nonuniformEXT(mat.emissiveTex)], fragUV).rgb;
    }

    outColor = vec4(color, albedo.a);
}
