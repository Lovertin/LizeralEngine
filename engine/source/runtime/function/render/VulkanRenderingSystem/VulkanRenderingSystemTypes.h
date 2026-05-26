#pragma once

#include "runtime/core/math/matrix4.h"
#include "runtime/core/math/vector2.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/vector4.h"
#include "runtime/function/render/rhi/vulkan/VulkanBuffer.h"

#include <cstdint>
#include <memory>

namespace Lizeral {

    struct GlobalFrameData {
        Matrix4x4 viewProj;
        Matrix4x4 invViewProj;
        Matrix4x4 prevViewProj;
        Vector3 cameraPos;
        float padding1;
        Vector3 lightDir;
        float lightIntensity;
        Vector3 lightColor;
        uint32_t frameIndex;
        Vector2 jitter;
        Vector2 padding2;
    };

    struct GPUPointLight {
        Vector4 posAndRadius;
        Vector4 colorAndIntensity;
    };

    struct VulkanFrameResource {
        std::unique_ptr<VulkanBuffer> buffer;
        uint64_t addr = 0;
        bool IsValid() const { return addr != 0; }
    };

    struct RTInstanceDesc {
        uint64_t vertexBuffer;
        uint64_t indexBuffer;
        uint64_t materialBuffer;
        uint64_t primitiveMaterialIdBuffer;
        uint32_t textureOffset;
        uint32_t materialCount;
        uint32_t padding[2];
    };

    struct GPUInstanceData {
        Matrix4x4 Model;
        Matrix4x4 prevModel;
    };

    struct MeshPassPushConstants {
        uint64_t frameDataAddr;
        uint64_t modelDataAddr;
        uint64_t vertexBuffer;
        uint64_t meshletBuffer;
        uint64_t indexBuffer;
        uint64_t boundsBuffer;
        uint64_t materialBuffer;
        uint32_t totalMeshlets;
        uint32_t textureOffset;
    };

    struct FullscreenPassPushConstants {
        uint64_t frameDataAddr;
        uint64_t instanceDescAddr;
        uint64_t pointLightsAddr;
        uint32_t stepSize;
        uint32_t pointLightCount;
    };

} // namespace Lizeral
