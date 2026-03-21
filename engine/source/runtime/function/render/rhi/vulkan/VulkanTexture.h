// source/runtime/function/render/rhi/vulkan/VulkanTexture.h
#pragma once
#include <vulkan/vulkan.h>
#include <vma/vk_mem_alloc.h>
#include "runtime/core/vendor/stb_image/stb_image.h"

#include <string>

namespace Lizeral {

    class VulkanDevice;
    class VulkanCommandPool;

    class VulkanTexture {
    public:
        VulkanTexture(VulkanDevice* device, VulkanCommandPool* commandPool, const std::string& filepath);

        VulkanTexture(VulkanDevice* device, VulkanCommandPool* commandPool, const unsigned char* buffer, int len);

        ~VulkanTexture();

        VulkanTexture(const VulkanTexture&) = delete;
        VulkanTexture& operator=(const VulkanTexture&) = delete;

        // RHI handle
        VkImage GetNativeImage() const { return m_image; }
        VkImageView GetImageView() const { return m_imageView; }
        VkSampler GetSampler() const { return m_sampler; }
        uint32_t GetMipLevels() const { return m_mipLevels; }

    private:
        VulkanDevice* m_device { nullptr };
        
        VkImage m_image { VK_NULL_HANDLE };
        VmaAllocation m_allocation { VK_NULL_HANDLE };
        VkImageView m_imageView { VK_NULL_HANDLE };
        VkSampler m_sampler { VK_NULL_HANDLE };

        uint32_t m_width { 0 };
        uint32_t m_height { 0 };
        uint32_t m_mipLevels { 1 };

        void createTextureImage(VulkanCommandPool* commandPool, const std::string& filepath);
        void createTextureImageFromPixels(VulkanCommandPool* commandPool,stbi_uc* pixels,int texWidth,int texHeight);
        void generateMipmaps(VkCommandBuffer cmd, VkImage image, VkFormat imageFormat, int32_t texWidth, int32_t texHeight, uint32_t mipLevels);
        void createImageView();
        void createTextureSampler();
    };

} // namespace Lizeral