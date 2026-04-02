#version 460
#extension GL_EXT_ray_query : require
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_scalar_block_layout : require
#extension GL_EXT_nonuniform_qualifier : require
#extension GL_EXT_shader_explicit_arithmetic_types_int64 : require
#extension GL_GOOGLE_include_directive : require

// Specialization Constants

layout(constant_id = 0) const int GI_QUALITY_LEVEL = 2; // 0: Off, 1: SSGI, 2: RTGI
layout(constant_id = 1) const int SHADOW_QUALITY = 1; // 0: hard shadow,1: soft shadow

#include "Common/SceneShared.glsl" 
#include "Common/Math.glsl"       
#include "Common/BRDF.glsl"    

// Shader IO & resource bine

layout(location = 0) in vec2 inUV;

layout(location = 0) out vec4 outDirectLight;
layout(location = 1) out vec4 outNoisyGI;

layout(binding = 0) uniform sampler2D samplerAlbedoMetallic;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;
layout(binding = 4) uniform accelerationStructureEXT topLevelAS;
layout(set = 1, binding = 0) uniform sampler2D bindlessTextures[];

#include "ssrm_lib.glsl"

vec3 TraceGlobalIlluminationRay(vec3 origin, vec3 direction, inout uint seed) {
    rayQueryEXT rayQuery;
    uint rayFlags = gl_RayFlagsOpaqueEXT;
    rayQueryInitializeEXT(rayQuery, topLevelAS, rayFlags, 0xFF, origin, 0.01, direction, 1000.0);
    while(rayQueryProceedEXT(rayQuery)) {}

    if (rayQueryGetIntersectionTypeEXT(rayQuery, true) == gl_RayQueryCommittedIntersectionTriangleEXT) {
        uint instanceID = rayQueryGetIntersectionInstanceCustomIndexEXT(rayQuery, true);
        uint primID = rayQueryGetIntersectionPrimitiveIndexEXT(rayQuery, true);
        vec2 bary = rayQueryGetIntersectionBarycentricsEXT(rayQuery, true);
        
        InstanceDescBuffer instances = InstanceDescBuffer(pc.instanceDescAddr);
        RTInstanceDesc desc = instances.instances[instanceID];
        IndexBuffer indices = IndexBuffer(desc.indexBuffer);
        uint i0 = indices.iData[primID * 3 + 0];
        uint i1 = indices.iData[primID * 3 + 1];
        uint i2 = indices.iData[primID * 3 + 2];
        VertexBuffer verts = VertexBuffer(desc.vertexBuffer);
        
        
        vec2 uv0 = vec2(verts.vData[i0 * 8 + 6], verts.vData[i0 * 8 + 7]);
        vec2 uv1 = vec2(verts.vData[i1 * 8 + 6], verts.vData[i1 * 8 + 7]);
        vec2 uv2 = vec2(verts.vData[i2 * 8 + 6], verts.vData[i2 * 8 + 7]);
        vec2 hitUV = uv0 * (1.0 - bary.x - bary.y) + uv1 * bary.x + uv2 * bary.y;
        vec3 hitAlbedo = vec3(1.0);
        float hitMetallic = 0.0;
        float hitRoughness = 1.0;
        if (desc.materialCount > 0u && desc.materialBuffer != 0ul && desc.primitiveMaterialIdBuffer != 0ul) {
            PrimitiveMaterialIdBuffer primMatIds = PrimitiveMaterialIdBuffer(desc.primitiveMaterialIdBuffer);
            MaterialBuffer materials = MaterialBuffer(desc.materialBuffer);
            uint materialId = min(primMatIds.materialId[primID], desc.materialCount - 1u);
            Material hitMaterial = materials.m[materialId];
            hitAlbedo = hitMaterial.baseColorFactor.rgb;
            hitMetallic = hitMaterial.metallicFactor;
            hitRoughness = hitMaterial.roughnessFactor;
            if (hitMaterial.albedoTex >= 0) {
                hitAlbedo *= texture(bindlessTextures[nonuniformEXT(hitMaterial.albedoTex)], hitUV).rgb;
            }
            if (hitMaterial.ormTex >= 0) {
                vec3 orm = texture(bindlessTextures[nonuniformEXT(hitMaterial.ormTex)], hitUV).rgb;
                hitRoughness *= orm.g;
                hitMetallic *= orm.b;
            }
        }
        hitMetallic = clamp(hitMetallic, 0.0, 1.0);
        hitRoughness = clamp(hitRoughness, 0.05, 1.0);
        
        vec3 bounceDiffuse = hitAlbedo * (1.0 - hitMetallic);
        float roughnessAtten = 1.0 - 0.5 * hitRoughness;
        vec3 bouncedLight = bounceDiffuse * pc.frameData.frame.lightColor * pc.frameData.frame.lightIntensity * 0.5 * roughnessAtten;
        return bouncedLight;
    }
    
    return mix(vec3(0.1, 0.2, 0.3), vec3(0.01, 0.02, 0.05), direction.y * 0.5 + 0.5); // Skybox Fallback
}

vec3 ComputeGeometricNormal(vec3 worldPos, vec3 shadingNormal) {
    vec3 dpdx = dFdx(worldPos);
    vec3 dpdy = dFdy(worldPos);
    vec3 geomNormal = normalize(cross(dpdx, dpdy));
    if (dot(geomNormal, shadingNormal) < 0.0) {
        geomNormal = -geomNormal;
    }
    return geomNormal;
}

vec3 StabilizeShadingNormal(vec3 shadingNormal, vec3 geomNormal) {
    vec3 n = normalize(shadingNormal);
    if (dot(n, geomNormal) < 0.0) {
        n = -n;
    }
    return n;
}

float ApplySpecularAA(vec3 normalWS, float roughness) {
    vec3 dndx = dFdx(normalWS);
    vec3 dndy = dFdy(normalWS);
    float variance = max(dot(dndx, dndx), dot(dndy, dndy));
    float aa = clamp(variance, 0.0, 1.0);
    return clamp(sqrt(roughness * roughness + aa), 0.05, 1.0);
}

void main() {
    float depth = texture(samplerDepth, inUV).r;
    if (depth <= 0.0) {
        vec3 bg = mix(vec3(0.1, 0.2, 0.3), vec3(0.01, 0.02, 0.05), inUV.y);
        outDirectLight = vec4(bg, 1.0);        
        outNoisyGI = vec4(0.0, 0.0, 0.0, 1.0); 
        return;
    }

    vec4 albedoMetallic = texture(samplerAlbedoMetallic, inUV);
    vec4 normalRoughness = texture(samplerNormalRoughness, inUV);
    vec3 albedo = albedoMetallic.rgb; 
    float metallic = clamp(albedoMetallic.a, 0.0, 0.9);
    vec3 normal = normalize(normalRoughness.rgb + vec3(0.0001));
    float roughness = clamp(normalRoughness.a, 0.05, 1.0); 

    // GBuffer depth is generated from jittered clip-space position in geometry pass.
    // Undo jitter in UV before reconstructing world position with unjittered invViewProj.
    vec2 unjitteredUV = inUV - pc.frameData.frame.jitter * 0.5;
    vec3 worldPos = ReconstructWorldPos(unjitteredUV, depth, pc.frameData.frame.invViewProj);
    vec3 viewDir = normalize(pc.frameData.frame.cameraPos - worldPos);
    vec3 geomNormal = ComputeGeometricNormal(worldPos, normal);

    uint baseSeed = uint(gl_FragCoord.x) * 1973u + uint(gl_FragCoord.y) * 9277u;
    uint giBaseSeed = baseSeed + pc.frameData.frame.frameIndex * 26699u;
    float rayBias = 0.002;
    vec3 rayOrigin = worldPos + geomNormal * rayBias;
    float tMin = rayBias;     
    uint rayFlags = gl_RayFlagsTerminateOnFirstHitEXT | gl_RayFlagsOpaqueEXT | gl_RayFlagsSkipClosestHitShaderEXT;
    uint shadowRayFlags = rayFlags | gl_RayFlagsCullBackFacingTrianglesEXT;

    vec3 finalDirectLight = vec3(0.0);
    vec3 mainLightDir = normalize(pc.frameData.frame.lightDir);
    
    float shadow = 1.0; 

    if (SHADOW_QUALITY == 1) {
        float lightRadius = 0.03; 
        vec3 lightUp = abs(mainLightDir.y) > 0.99 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
        vec3 lightRight = normalize(cross(lightUp, mainLightDir));
        lightUp = cross(mainLightDir, lightRight);

        int NUM_SAMPLES = 8;
        float shadowSum = 0.0; 
        float baseAngle = rand(baseSeed) * 2.0 * 3.14159265;
        
        for(int i = 0; i < NUM_SAMPLES; i++) {
            float r = sqrt((float(i) + 0.5) / float(NUM_SAMPLES)) * lightRadius;
            float theta = float(i) * 2.39996323 + baseAngle;
            vec2 diskOffset = vec2(r * cos(theta), r * sin(theta));
            vec3 jitteredRayDir = normalize(mainLightDir + lightRight * diskOffset.x + lightUp * diskOffset.y);

            float ndotlShadow = abs(dot(geomNormal, jitteredRayDir));
            float shadowBias = mix(0.004, 0.0005, ndotlShadow);
            vec3 shadowOrigin = worldPos + geomNormal * shadowBias + jitteredRayDir * (shadowBias * 0.2);

            rayQueryEXT rayQuery;
            rayQueryInitializeEXT(rayQuery, topLevelAS, shadowRayFlags, 0xFF, shadowOrigin, shadowBias, jitteredRayDir, 1000.0);
            while(rayQueryProceedEXT(rayQuery)) {}
            if (rayQueryGetIntersectionTypeEXT(rayQuery, true) == gl_RayQueryCommittedIntersectionNoneEXT) {
                shadowSum += 1.0;
            }
        }
        shadow = shadowSum / float(NUM_SAMPLES);
    } else if(SHADOW_QUALITY == 0){
        float ndotlShadow = abs(dot(geomNormal, mainLightDir));
        float shadowBias = mix(0.004, 0.0005, ndotlShadow);
        vec3 shadowOrigin = worldPos + geomNormal * shadowBias + mainLightDir * (shadowBias * 0.2);

        rayQueryEXT rayQuery;
        rayQueryInitializeEXT(rayQuery, topLevelAS, shadowRayFlags, 0xFF, shadowOrigin, shadowBias, mainLightDir, 1000.0);
        while(rayQueryProceedEXT(rayQuery)) {}
        
        if (rayQueryGetIntersectionTypeEXT(rayQuery, true) != gl_RayQueryCommittedIntersectionNoneEXT) {
            shadow = 0.0; 
        } else {
            shadow = 1.0; 
        }
    }

    float NdotL_main;
    vec3 mainBrdf = CalculateBRDF(normal, viewDir, mainLightDir, albedo, metallic, roughness, NdotL_main);
    finalDirectLight += mainBrdf * pc.frameData.frame.lightColor * pc.frameData.frame.lightIntensity * NdotL_main * shadow;

    if (pc.pointLightsAddr != 0 && pc.pointLightCount > 0) {
        PointLightBuffer plBuf = PointLightBuffer(pc.pointLightsAddr);
        
        for (uint i = 0; i < pc.pointLightCount; i++) {
            GPUPointLight pl = plBuf.lights[i];
            vec3 lightPos = pl.posAndRadius.xyz;
            float radius = pl.posAndRadius.w;
            vec3 lightColor = pl.colorAndIntensity.rgb;
            float intensity = pl.colorAndIntensity.a;

            vec3 L = lightPos - worldPos;
            float distance = length(L);

            if (distance < radius) {
                L /= distance; 

                float falloff = CalculatePointLightFalloff(distance, radius);

                const bool ENABLE_POINT_LIGHT_SHADOWS = false;
                float plShadow = 1.0;
                if (ENABLE_POINT_LIGHT_SHADOWS) {
                    float ndotlShadow = abs(dot(geomNormal, L));
                    float shadowBias = mix(0.004, 0.0005, ndotlShadow);
                    vec3 shadowOrigin = worldPos + geomNormal * shadowBias + L * (shadowBias * 0.2);

                    rayQueryEXT rq;
                    rayQueryInitializeEXT(rq, topLevelAS, shadowRayFlags, 0xFF, shadowOrigin, shadowBias, L, distance);
                    while(rayQueryProceedEXT(rq)) {}
                    if (rayQueryGetIntersectionTypeEXT(rq, true) != gl_RayQueryCommittedIntersectionNoneEXT) {
                        plShadow = 0.0;
                    }
                }

                float NdotL_pl;
                vec3 plBrdf = CalculateBRDF(normal, viewDir, L, albedo, metallic, roughness, NdotL_pl);
                
                finalDirectLight += plBrdf * lightColor * intensity * NdotL_pl * falloff * plShadow;
            }
        }
    }

    outDirectLight = vec4(finalDirectLight, 1.0);

    // Global Illumination multi level
    vec3 indirectLight = vec3(0.0);

    vec3 w_tbn = normal;
    vec3 u_tbn = normalize(cross((abs(w_tbn.x) > 0.1 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0)), w_tbn));
    vec3 v_tbn = cross(w_tbn, u_tbn);
    mat3 tbn_hemisphere = mat3(u_tbn, v_tbn, w_tbn);

    int numGisamples = 1; 
    vec3 giSum = vec3(0.0);
    uint giSeed = giBaseSeed + pc.frameData.frame.frameIndex * 161803u; 

    if (GI_QUALITY_LEVEL == 2) {
        // 2 level RTGI
        for(int i = 0; i < numGisamples; i++) {
            float r1 = rand(giSeed);
            float r2 = rand(giSeed);
            float r = sqrt(r1);
            float phi = 2.0 * 3.14159265 * r2;
            vec3 localDir = vec3(r * cos(phi), r * sin(phi), sqrt(max(0.0, 1.0 - r1)));
            vec3 giDirection = normalize(tbn_hemisphere * localDir);
            
            vec3 col = TraceGlobalIlluminationRay(rayOrigin, giDirection, giSeed);
            giSum += min(col, vec3(2.0)); // avoid Fireflies
        }
        
        vec3 bouncedColor = (giSum / float(numGisamples)) * 3.14159265;
        float giMultiplier = sqrt(pc.frameData.frame.lightIntensity) * 0.02;
        indirectLight = bouncedColor * giMultiplier;

    } else if (GI_QUALITY_LEVEL == 1) {
        float r1 = rand(giSeed);
        float r2 = rand(giSeed);
        float r = sqrt(r1);
        float phi = 2.0 * 3.14159265 * r2;
        vec3 localDir = vec3(r * cos(phi), r * sin(phi), sqrt(max(0.0, 1.0 - r1)));
        vec3 giDirection = normalize(tbn_hemisphere * localDir);

        float hitMask = 0.0;
        vec3 hitColor = ScreenSpaceRayMarch(worldPos + geomNormal * rayBias, giDirection, hitMask);

        if (hitMask > 0.5) {
            indirectLight = hitColor * albedo * pc.frameData.frame.lightColor * sqrt(pc.frameData.frame.lightIntensity) * 0.05;
        } else {
            indirectLight = albedo * 0.02 * pc.frameData.frame.lightColor;
        }

    } else {
        indirectLight = albedo * 0.01;
    }
    
    if (GI_QUALITY_LEVEL > 0) {
        outNoisyGI = vec4(indirectLight, 1.0);
    } else {
        outDirectLight += vec4(indirectLight, 0.0);
        outNoisyGI = vec4(0.0, 0.0, 0.0, 1.0);
    }

    outNoisyGI = vec4(indirectLight, 1.0);
}


