#pragma once

#include <vma/vk_mem_alloc.h>
#include <vulkan/vulkan.h>

#include <cstdint>

namespace Lizeral {

    class VulkanDevice;

    struct GBufferAttachment {
        VkImage image = VK_NULL_HANDLE;
        VmaAllocation allocation = VK_NULL_HANDLE;
        VkImageView view = VK_NULL_HANDLE;
        VkFormat format = VK_FORMAT_UNDEFINED;
    };

    class RenderFrameTargets {
    public:
        void Initialize(VulkanDevice* device, uint32_t width, uint32_t height);
        void Shutdown();
        void Resize(uint32_t width, uint32_t height);

        uint32_t GetWidth() const { return m_width; }
        uint32_t GetHeight() const { return m_height; }
        VkSampler GetSampler() const { return m_sampler; }

        GBufferAttachment albedoMetallic;
        GBufferAttachment normalRoughness;
        GBufferAttachment depth;
        GBufferAttachment velocity;
        GBufferAttachment directLight;
        GBufferAttachment noisyGI;
        GBufferAttachment denoisedGI;
        GBufferAttachment denoisedGITemp;
        GBufferAttachment giHistory[2];
        GBufferAttachment momentsHistory[2];
        GBufferAttachment history[2];

    private:
        GBufferAttachment CreateAttachment(VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect);
        void CreateAttachments();
        void DestroyAttachment(GBufferAttachment& attachment);
        void DestroyAttachments();
        void CreateSampler();
        void DestroySampler();

        VulkanDevice* m_device { nullptr };
        uint32_t m_width = 0;
        uint32_t m_height = 0;
        VkSampler m_sampler { VK_NULL_HANDLE };
    };

} // namespace Lizeral
