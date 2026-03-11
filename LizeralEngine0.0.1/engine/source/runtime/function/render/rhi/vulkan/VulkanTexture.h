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
        // 构造函数：负责从硬盘读取图片、计算 Mipmap 层级、并传送到 GPU 显存
        VulkanTexture(VulkanDevice* device, VulkanCommandPool* commandPool, const std::string& filepath);

        VulkanTexture(VulkanDevice* device, VulkanCommandPool* commandPool, const unsigned char* buffer, int len);

        ~VulkanTexture();

        // 禁用拷贝
        VulkanTexture(const VulkanTexture&) = delete;
        VulkanTexture& operator=(const VulkanTexture&) = delete;

        // --- 核心获取接口 ---
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

        // 内部流程辅助函数
        void createTextureImage(VulkanCommandPool* commandPool, const std::string& filepath);
        void createTextureImageFromPixels(VulkanCommandPool* commandPool,stbi_uc* pixels,int texWidth,int texHeight);
        void generateMipmaps(VkCommandBuffer cmd, VkImage image, VkFormat imageFormat, int32_t texWidth, int32_t texHeight, uint32_t mipLevels);
        void createImageView();
        void createTextureSampler();
    };

} // namespace Lizeral