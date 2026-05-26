#pragma once

#include "runtime/core/math/matrix4.h"
#include "runtime/function/ecs/entity.h"
#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystemTypes.h"
#include "runtime/function/render/rhi/vulkan/VulkanBuffer.h"
#include "runtime/function/Vulkan_res_type/VulkanModelResource.h"

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>
#include <vulkan/vulkan.h>

namespace Lizeral {

    class RenderResourceCache;
    class RenderSceneSnapshot;

    struct RenderMeshDrawItem {
        Entity entity = null_entity;
        const VulkanModelResource* resource = nullptr;
        Matrix4x4 modelMatrix;
        uint64_t materialBuffer = 0;
        uint32_t instanceIndex = 0;
    };

    struct RenderTransparentDrawItem {
        uint32_t meshDrawIndex = 0;
        float distanceSq = 0.0f;
    };

    struct RenderDrawPacket {
        std::vector<RenderMeshDrawItem> meshDraws;
        std::vector<RenderTransparentDrawItem> transparentDraws;
        std::vector<VkAccelerationStructureInstanceKHR> tlasInstances;
        std::vector<RTInstanceDesc> rtInstanceDescs;
        std::unordered_map<uint32_t, uint32_t> entityToInstanceIndex;
    };

    RenderDrawPacket BuildRenderDrawPacket(const RenderSceneSnapshot& scene, RenderResourceCache& resourceCache);

    void WriteRenderDrawPacketInstanceData(
        const RenderDrawPacket& packet,
        bool isFirstFrame,
        std::unordered_map<uint32_t, Matrix4x4>& previousModelMatrices,
        VulkanBuffer& globalInstanceBuffer
    );

} // namespace Lizeral
