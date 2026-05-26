#include "RenderFrameTargets.h"

#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"

namespace Lizeral {

    void RenderFrameTargets::Initialize(VulkanDevice* device, uint32_t width, uint32_t height) {
        m_device = device;
        m_width = width;
        m_height = height;
        CreateAttachments();
        CreateSampler();
    }

    void RenderFrameTargets::Shutdown() {
        DestroyAttachments();
        DestroySampler();
        m_device = nullptr;
        m_width = 0;
        m_height = 0;
    }

    void RenderFrameTargets::Resize(uint32_t width, uint32_t height) {
        if (width == 0 || height == 0) {
            return;
        }

        if (width == m_width && height == m_height) {
            return;
        }

        DestroyAttachments();
        m_width = width;
        m_height = height;
        CreateAttachments();
    }

    GBufferAttachment RenderFrameTargets::CreateAttachment(VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect) {
        GBufferAttachment attachment;
        attachment.format = format;

        VkImageCreateInfo imageInfo{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
        imageInfo.imageType = VK_IMAGE_TYPE_2D;
        imageInfo.extent = {m_width, m_height, 1};
        imageInfo.mipLevels = 1;
        imageInfo.arrayLayers = 1;
        imageInfo.format = format;
        imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
        imageInfo.usage = usage | VK_IMAGE_USAGE_SAMPLED_BIT;
        imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
        imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

        VmaAllocationCreateInfo allocInfo{};
        allocInfo.usage = VMA_MEMORY_USAGE_GPU_ONLY;

        vmaCreateImage(m_device->GetAllocator(), &imageInfo, &allocInfo, &attachment.image, &attachment.allocation, nullptr);

        VkImageViewCreateInfo viewInfo{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
        viewInfo.image = attachment.image;
        viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        viewInfo.format = format;
        viewInfo.subresourceRange.aspectMask = aspect;
        viewInfo.subresourceRange.levelCount = 1;
        viewInfo.subresourceRange.layerCount = 1;

        vkCreateImageView(m_device->GetNativeDevice(), &viewInfo, nullptr, &attachment.view);
        return attachment;
    }

    void RenderFrameTargets::CreateAttachments() {
        albedoMetallic = CreateAttachment(VK_FORMAT_R8G8B8A8_UNORM, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        normalRoughness = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        depth = CreateAttachment(VK_FORMAT_D32_SFLOAT, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT, VK_IMAGE_ASPECT_DEPTH_BIT);
        velocity = CreateAttachment(VK_FORMAT_R16G16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

        directLight = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        noisyGI = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        denoisedGI = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        denoisedGITemp = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

        giHistory[0] = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        giHistory[1] = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

        momentsHistory[0] = CreateAttachment(VK_FORMAT_R16G16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        momentsHistory[1] = CreateAttachment(VK_FORMAT_R16G16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

        history[0] = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        history[1] = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
    }

    void RenderFrameTargets::DestroyAttachment(GBufferAttachment& attachment) {
        if (m_device == nullptr) {
            return;
        }

        if (attachment.view != VK_NULL_HANDLE) {
            vkDestroyImageView(m_device->GetNativeDevice(), attachment.view, nullptr);
            attachment.view = VK_NULL_HANDLE;
        }
        if (attachment.image != VK_NULL_HANDLE && attachment.allocation != VK_NULL_HANDLE) {
            vmaDestroyImage(m_device->GetAllocator(), attachment.image, attachment.allocation);
            attachment.image = VK_NULL_HANDLE;
            attachment.allocation = VK_NULL_HANDLE;
        }
        attachment.format = VK_FORMAT_UNDEFINED;
    }

    void RenderFrameTargets::DestroyAttachments() {
        DestroyAttachment(albedoMetallic);
        DestroyAttachment(normalRoughness);
        DestroyAttachment(depth);
        DestroyAttachment(velocity);
        DestroyAttachment(directLight);
        DestroyAttachment(noisyGI);
        DestroyAttachment(denoisedGI);
        DestroyAttachment(denoisedGITemp);
        DestroyAttachment(giHistory[0]);
        DestroyAttachment(giHistory[1]);
        DestroyAttachment(momentsHistory[0]);
        DestroyAttachment(momentsHistory[1]);
        DestroyAttachment(history[0]);
        DestroyAttachment(history[1]);
    }

    void RenderFrameTargets::CreateSampler() {
        VkSamplerCreateInfo samplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
        samplerInfo.magFilter = VK_FILTER_LINEAR;
        samplerInfo.minFilter = VK_FILTER_LINEAR;
        samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
        samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
        samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
        samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
        vkCreateSampler(m_device->GetNativeDevice(), &samplerInfo, nullptr, &m_sampler);
    }

    void RenderFrameTargets::DestroySampler() {
        if (m_device != nullptr && m_sampler != VK_NULL_HANDLE) {
            vkDestroySampler(m_device->GetNativeDevice(), m_sampler, nullptr);
            m_sampler = VK_NULL_HANDLE;
        }
    }

} // namespace Lizeral
