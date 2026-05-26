#pragma once

#include "runtime/function/ecs/entity.h"
#include "runtime/function/render/rhi/vulkan/VulkanBuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanTexture.h"
#include "runtime/function/Vulkan_res_type/VulkanModelResource.h"

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <vulkan/vulkan.h>

namespace Lizeral {

    class VulkanCommandPool;
    class VulkanDevice;
    class VulkanModelComponent;
    struct VulkanMaterialSlotOverride;

    class RenderResourceCache {
    public:
        static constexpr uint32_t kMaxBindlessTextures = 1024;

        void Initialize(VulkanDevice* device, VulkanCommandPool* resourceCommandPool);
        void SetBindlessDescriptorSet(VkDescriptorSet descriptorSet);
        void Clear();

        VulkanModelResource& GetOrLoadModel(const std::string& path);
        uint32_t GetOrLoadTextureIndex(const std::string& texturePath);
        uint64_t ResolveMaterialBufferAddress(Entity entity, const VulkanModelComponent& modelComp, const VulkanModelResource& modelResource);
        uint64_t ResolveMaterialBufferAddress(
            Entity entity,
            const std::string& modelAssetPath,
            uint32_t resourceRevision,
            const std::vector<VulkanMaterialSlotOverride>& materialOverrides,
            const VulkanModelResource& modelResource
        );

        const std::vector<VkDescriptorImageInfo>& GetGlobalImageInfos() const { return m_globalImageInfos; }

    private:
        struct OverrideMaterialBufferCacheEntry {
            std::string modelAssetPath;
            uint32_t componentRevision = 0;
            std::unique_ptr<VulkanBuffer> materialBuffer;
            uint64_t materialBufferAddress = 0;
        };

        template <typename T>
        std::unique_ptr<VulkanBuffer> CreateBDABuffer(const std::vector<T>& data) {
            if (data.empty()) {
                return nullptr;
            }

            VkDeviceSize bufferSize = data.size() * sizeof(T);
            auto buffer = std::make_unique<VulkanBuffer>(
                m_device,
                bufferSize,
                VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT | VK_BUFFER_USAGE_ACCELERATION_STRUCTURE_BUILD_INPUT_READ_ONLY_BIT_KHR,
                VMA_MEMORY_USAGE_CPU_TO_GPU
            );

            buffer->WriteData(data.data(), bufferSize);
            return buffer;
        }

        void ApplyMaterialOverride(MaterialData& targetMaterial, const VulkanMaterialSlotOverride& overrideData);
        void UpdateBindlessTextureDescriptor(uint32_t textureIndex);

        VulkanDevice* m_device { nullptr };
        VulkanCommandPool* m_resourceCommandPool { nullptr };
        VkDescriptorSet m_bindlessDescriptorSet { VK_NULL_HANDLE };

        std::unordered_map<std::string, VulkanModelResource> m_modelCache;
        std::vector<std::unique_ptr<VulkanTexture>> m_globalTextures;
        std::vector<VkDescriptorImageInfo> m_globalImageInfos;
        std::unordered_map<std::string, uint32_t> m_globalTexturePathCache;
        std::unordered_map<uint32_t, OverrideMaterialBufferCacheEntry> m_overrideMaterialBuffers;
    };

} // namespace Lizeral
