// shader/Common/SceneShared.glsl
#ifndef SCENE_SHARED_GLSL
#define SCENE_SHARED_GLSL

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

struct GPUPointLight {
    vec4 posAndRadius;     
    vec4 colorAndIntensity; 
};

layout(buffer_reference, scalar, buffer_reference_align = 16) readonly buffer PointLightBuffer { 
    GPUPointLight lights[]; 
};

struct RTInstanceDesc {
    uint64_t vertexBuffer;
    uint64_t indexBuffer;
    uint64_t materialBuffer;
    uint textureOffset;
    uint pad0, pad1, pad2;
};

layout(push_constant) uniform PushConstants {
    FrameDataBuffer frameData;
    uint64_t instanceDescAddr;       
    uint64_t pointLightsAddr;  
    uint stepSize;             
    uint pointLightCount;    
} pc;

// BDA 
layout(buffer_reference, scalar) readonly buffer VertexBuffer { float vData[]; };
layout(buffer_reference, scalar) readonly buffer IndexBuffer { uint iData[]; };
layout(buffer_reference, scalar) readonly buffer InstanceDescBuffer { RTInstanceDesc instances[]; };

#endif // SCENE_SHARED_GLSL