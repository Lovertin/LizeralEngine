// source/runtime/function/render/rhi/vulkan/VulkanTexture.cpp
#include "VulkanTexture.h"
#include "VulkanDevice.h"
#include "VulkanCommandPool.h"
#include "VulkanCommandBuffer.h"
#include "VulkanBuffer.h"

// 引入 stb_image 读取图片
// #define STB_IMAGE_IMPLEMENTATION

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace Lizeral {

    VulkanTexture::VulkanTexture(VulkanDevice* device, VulkanCommandPool* commandPool, const std::string& filepath)
        : m_device(device) {
        
        createTextureImage(commandPool, filepath);
        createImageView();
        createTextureSampler();
    }

    // ★ 新增的内存加载构造函数
    VulkanTexture::VulkanTexture(VulkanDevice* device, VulkanCommandPool* commandPool, const unsigned char* buffer, int len)
        : m_device(device) {
        
        int texWidth, texHeight, texChannels;
        
        // 尝试解析图片
        stbi_uc* pixels = stbi_load_from_memory(buffer, len, &texWidth, &texHeight, &texChannels, STBI_rgb_alpha);
        
        if (!pixels) {
            // ★ 核心修复：不要崩溃！打印出 stb_image 真实的失败原因
            std::cerr << "[VulkanTexture] WARNING: Failed to decode image! Reason: " << stbi_failure_reason() << std::endl;
            std::cerr << "[VulkanTexture] Creating 1x1 fallback white texture instead." << std::endl;
            
            // ★ 手工捏一张 1x1 的纯白贴图救场！
            texWidth = 1;
            texHeight = 1;
            // malloc 分配 4 个字节，因为后面会统一调用 stbi_image_free(它底层也是调用 free)
            pixels = (stbi_uc*)malloc(4); 
            pixels[0] = 255; // R
            pixels[1] = 255; // G
            pixels[2] = 255; // B
            pixels[3] = 255; // A
        } else {
            std::cout << "[VulkanTexture] Loaded embedded image from GLB (" << texWidth << "x" << texHeight << ")" << std::endl;
        }

        // 继续走后续的显卡提交流程
        createTextureImageFromPixels(commandPool, pixels, texWidth, texHeight);
        createImageView();
        createTextureSampler();
    }

    VulkanTexture::~VulkanTexture() {
        VkDevice device = m_device->GetNativeDevice();
        
        if (m_sampler != VK_NULL_HANDLE) {
            vkDestroySampler(device, m_sampler, nullptr);
        }
        if (m_imageView != VK_NULL_HANDLE) {
            vkDestroyImageView(device, m_imageView, nullptr);
        }
        if (m_image != VK_NULL_HANDLE && m_allocation != VK_NULL_HANDLE) {
            vmaDestroyImage(m_device->GetAllocator(), m_image, m_allocation);
        }
    }

    void VulkanTexture::createTextureImage(VulkanCommandPool* commandPool, const std::string& filepath) {
        // 1. 使用 stb_image 读取像素数据 (强制要求 4 通道 RGBA)
        int texWidth, texHeight, texChannels;
        stbi_uc* pixels = stbi_load(filepath.c_str(), &texWidth, &texHeight, &texChannels, STBI_rgb_alpha);
        
        if (!pixels) {
            throw std::runtime_error("Failed to load texture image: " + filepath);
        }
        std::cout << "[VulkanTexture] Loaded image: " << filepath << " (" << texWidth << "x" << texHeight << ")" << std::endl;

        m_width = static_cast<uint32_t>(texWidth);
        m_height = static_cast<uint32_t>(texHeight);
        VkDeviceSize imageSize = m_width * m_height * 4; // 4 bytes per pixel (RGBA)

        // 计算这图能生成多少层 Mipmap（公式：log2(max(宽, 高)) + 1）
        m_mipLevels = static_cast<uint32_t>(std::floor(std::log2(std::max(m_width, m_height)))) + 1;

        // 2. 创建 Staging Buffer (中转站，位于 CPU 和 GPU 之间)
        VulkanBuffer stagingBuffer(m_device, imageSize, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, VMA_MEMORY_USAGE_CPU_ONLY);
        stagingBuffer.WriteData(pixels, imageSize); // 把像素拷贝进去
        stbi_image_free(pixels); // 释放主内存里的像素

        // 3. 通过 VMA 创建真正的 GPU 显存贴图
        VkImageCreateInfo imageInfo{};
        imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
        imageInfo.imageType = VK_IMAGE_TYPE_2D;
        imageInfo.extent.width = m_width;
        imageInfo.extent.height = m_height;
        imageInfo.extent.depth = 1;
        imageInfo.mipLevels = m_mipLevels;
        imageInfo.arrayLayers = 1;
        imageInfo.format = VK_FORMAT_R8G8B8A8_SRGB; // Albedo 贴图必须是 SRGB 空间！
        imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
        imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        // 用途：作为拷贝目标、作为 Mipmap 的拷贝源、作为 Shader 采样的图
        imageInfo.usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
        imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;

        VmaAllocationCreateInfo allocInfo{};
        allocInfo.usage = VMA_MEMORY_USAGE_GPU_ONLY;

        if (vmaCreateImage(m_device->GetAllocator(), &imageInfo, &allocInfo, &m_image, &m_allocation, nullptr) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create image via VMA!");
        }

        // ========================================================
        // 4. 录制一次性指令，完成格式转换、数据搬运和 Mipmap 生成
        // ========================================================
        VulkanCommandBuffer cmdBuf(m_device, commandPool);
        cmdBuf.Begin(VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT);
        VkCommandBuffer cmd = cmdBuf.GetNativeBuffer();

        // (1) 把新图片的布局从 UNDEFINED 转换为 "准备接收拷贝 (TRANSFER_DST)"
        VkImageMemoryBarrier barrier{};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = m_image;
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = m_mipLevels; // 改变所有层的布局
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;
        barrier.srcAccessMask = 0;
        barrier.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;

        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0,
            0, nullptr, 0, nullptr, 1, &barrier);

        // (2) 把 Staging Buffer 里的像素拷贝到贴图的第 0 层
        VkBufferImageCopy region{};
        region.bufferOffset = 0;
        region.bufferRowLength = 0;
        region.bufferImageHeight = 0;
        region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        region.imageSubresource.mipLevel = 0; // 只拷贝到最高清的第 0 层
        region.imageSubresource.baseArrayLayer = 0;
        region.imageSubresource.layerCount = 1;
        region.imageOffset = {0, 0, 0};
        region.imageExtent = {m_width, m_height, 1};

        vkCmdCopyBufferToImage(cmd, stagingBuffer.GetNativeBuffer(), m_image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &region);

        // (3) 生成 Mipmap！(在这个函数里，图片的状态会顺便被转成 SHADER_READ_ONLY)
        generateMipmaps(cmd, m_image, VK_FORMAT_R8G8B8A8_SRGB, texWidth, texHeight, m_mipLevels);

        // 提交并暴力等待完成 (因为离开这个函数 stagingBuffer 就会被销毁)
        cmdBuf.End();
        cmdBuf.SubmitAndIdle();
    }

    void VulkanTexture::createTextureImageFromPixels(VulkanCommandPool* commandPool, stbi_uc* pixels, int texWidth, int texHeight) {
        
        // ★ 核心修复：必须第一时间赋值给成员变量！
        m_width = static_cast<uint32_t>(texWidth);
        m_height = static_cast<uint32_t>(texHeight);
        
        VkDeviceSize imageSize = m_width * m_height * 4; // 4 bytes per pixel (RGBA)

        // 计算这图能生成多少层 Mipmap（公式：log2(max(宽, 高)) + 1）
        m_mipLevels = static_cast<uint32_t>(std::floor(std::log2(std::max(m_width, m_height)))) + 1;

        // 2. 创建 Staging Buffer (中转站，位于 CPU 和 GPU 之间)
        VulkanBuffer stagingBuffer(m_device, imageSize, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, VMA_MEMORY_USAGE_CPU_ONLY);
        stagingBuffer.WriteData(pixels, imageSize);

        stbi_image_free(pixels); 

        // 3. 通过 VMA 创建真正的 GPU 显存贴图
        VkImageCreateInfo imageInfo{};
        imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
        imageInfo.imageType = VK_IMAGE_TYPE_2D;
        imageInfo.extent.width = m_width;   // 使用成员变量更统一
        imageInfo.extent.height = m_height;
        imageInfo.extent.depth = 1;
        imageInfo.mipLevels = m_mipLevels;
        imageInfo.arrayLayers = 1;
        imageInfo.format = VK_FORMAT_R8G8B8A8_SRGB; // Albedo 贴图必须是 SRGB 空间！
        imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
        imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

        imageInfo.usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
        imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;

        VmaAllocationCreateInfo allocInfo{};
        allocInfo.usage = VMA_MEMORY_USAGE_GPU_ONLY;

        if (vmaCreateImage(m_device->GetAllocator(), &imageInfo, &allocInfo, &m_image, &m_allocation, nullptr) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create image via VMA!");
        }

        VulkanCommandBuffer cmdBuf(m_device, commandPool);
        cmdBuf.Begin(VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT);
        VkCommandBuffer cmd = cmdBuf.GetNativeBuffer();

        VkImageMemoryBarrier barrier{};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = m_image;
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = m_mipLevels; 
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;
        barrier.srcAccessMask = 0;
        barrier.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;

        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0,
            0, nullptr, 0, nullptr, 1, &barrier);

        VkBufferImageCopy region{};
        region.bufferOffset = 0;
        region.bufferRowLength = 0;
        region.bufferImageHeight = 0;
        region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        region.imageSubresource.mipLevel = 0; 
        region.imageSubresource.baseArrayLayer = 0;
        region.imageSubresource.layerCount = 1;
        region.imageOffset = {0, 0, 0};
        region.imageExtent = {m_width, m_height, 1}; 

        vkCmdCopyBufferToImage(cmd, stagingBuffer.GetNativeBuffer(), m_image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &region);

        generateMipmaps(cmd, m_image, VK_FORMAT_R8G8B8A8_SRGB, m_width, m_height, m_mipLevels);

        cmdBuf.End();
        cmdBuf.SubmitAndIdle();
    }

    void VulkanTexture::generateMipmaps(VkCommandBuffer cmd, VkImage image, VkFormat imageFormat, int32_t texWidth, int32_t texHeight, uint32_t mipLevels) {
        VkImageMemoryBarrier barrier{};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.image = image;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;
        barrier.subresourceRange.levelCount = 1;

        int32_t mipWidth = texWidth;
        int32_t mipHeight = texHeight;

        // 循环把上一层 (i-1) 缩小一半拷贝到当前层 (i)
        for (uint32_t i = 1; i < mipLevels; i++) {
            
            // 1. 把上一层 (i-1) 从 DST 变成 SRC (准备被读取)
            barrier.subresourceRange.baseMipLevel = i - 1;
            barrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
            barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
            barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
            barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
            vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0,
                0, nullptr, 0, nullptr, 1, &barrier);

            // 2. 发起 Blit (缩放拷贝) 命令
            VkImageBlit blit{};
            blit.srcOffsets[0] = { 0, 0, 0 };
            blit.srcOffsets[1] = { mipWidth, mipHeight, 1 };
            blit.srcSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
            blit.srcSubresource.mipLevel = i - 1;
            blit.srcSubresource.baseArrayLayer = 0;
            blit.srcSubresource.layerCount = 1;
            
            blit.dstOffsets[0] = { 0, 0, 0 };
            blit.dstOffsets[1] = { mipWidth > 1 ? mipWidth / 2 : 1, mipHeight > 1 ? mipHeight / 2 : 1, 1 };
            blit.dstSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
            blit.dstSubresource.mipLevel = i;
            blit.dstSubresource.baseArrayLayer = 0;
            blit.dstSubresource.layerCount = 1;

            vkCmdBlitImage(cmd,
                image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                1, &blit, VK_FILTER_LINEAR); // 线性插值

            // 3. 上一层 (i-1) 用完了，直接把它转成 Shader 只读状态 (SHADER_READ_ONLY)
            barrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
            barrier.newLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
            barrier.srcAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
            barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
            vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0,
                0, nullptr, 0, nullptr, 1, &barrier);

            if (mipWidth > 1) mipWidth /= 2;
            if (mipHeight > 1) mipHeight /= 2;
        }

        // 把最后一层也转成 Shader 只读状态
        barrier.subresourceRange.baseMipLevel = mipLevels - 1;
        barrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
        barrier.newLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0,
            0, nullptr, 0, nullptr, 1, &barrier);
    }

    void VulkanTexture::createImageView() {
        VkImageViewCreateInfo viewInfo{};
        viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        viewInfo.image = m_image;
        viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        viewInfo.format = VK_FORMAT_R8G8B8A8_SRGB; // 保持一致
        viewInfo.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        viewInfo.subresourceRange.baseMipLevel = 0;
        viewInfo.subresourceRange.levelCount = m_mipLevels; // 相框能看到所有的层级
        viewInfo.subresourceRange.baseArrayLayer = 0;
        viewInfo.subresourceRange.layerCount = 1;

        if (vkCreateImageView(m_device->GetNativeDevice(), &viewInfo, nullptr, &m_imageView) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create texture image view!");
        }
    }

    void VulkanTexture::createTextureSampler() {
        VkSamplerCreateInfo samplerInfo{};
        samplerInfo.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
        samplerInfo.magFilter = VK_FILTER_LINEAR; // 放大时线性插值
        samplerInfo.minFilter = VK_FILTER_LINEAR; // 缩小时线性插值
        samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_REPEAT; // 纹理坐标超出时重复
        samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_REPEAT;
        samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_REPEAT;

        // 开启各向异性过滤 (让地面这种远处的纹理更清晰)
        samplerInfo.anisotropyEnable = VK_TRUE; 
        samplerInfo.maxAnisotropy = 16.0f; // 显卡一般最高支持 16X

        samplerInfo.borderColor = VK_BORDER_COLOR_INT_OPAQUE_BLACK;
        samplerInfo.unnormalizedCoordinates = VK_FALSE;
        samplerInfo.compareEnable = VK_FALSE;
        samplerInfo.compareOp = VK_COMPARE_OP_ALWAYS;

        // Mipmap 采样设置
        samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_LINEAR;
        samplerInfo.mipLodBias = 0.0f;
        samplerInfo.minLod = 0.0f;
        samplerInfo.maxLod = static_cast<float>(m_mipLevels);

        if (vkCreateSampler(m_device->GetNativeDevice(), &samplerInfo, nullptr, &m_sampler) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create texture sampler!");
        }
    }

} // namespace Lizeral