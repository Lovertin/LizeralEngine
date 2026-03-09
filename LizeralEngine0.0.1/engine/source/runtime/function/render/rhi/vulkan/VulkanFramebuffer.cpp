#include "VulkanFramebuffer.h"
#include "VulkanDevice.h"
#include "VulkanRenderPass.h"
#include "VulkanSwapchain.h"
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanFramebuffer::VulkanFramebuffer(VulkanDevice* device, VulkanRenderPass* renderPass, VulkanSwapchain* swapchain) 
        : m_device(device) {
        
        const auto& imageViews = swapchain->GetImageViews();
        VkExtent2D extent = swapchain->GetExtent();
        
        m_framebuffers.resize(imageViews.size());

        // 为 Swapchain 里的每一张 ImageView 创建一个对应的 Framebuffer
        for (size_t i = 0; i < imageViews.size(); i++) {
            VkImageView attachments[] = {
                imageViews[i] // 靶子：当前正在处理的这张 Swapchain 图像
            };

            VkFramebufferCreateInfo framebufferInfo{};
            framebufferInfo.sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO;
            framebufferInfo.renderPass = renderPass->GetNativeRenderPass(); // 绑定说明书
            framebufferInfo.attachmentCount = 1;
            framebufferInfo.pAttachments = attachments; // 绑定靶子
            framebufferInfo.width = extent.width;
            framebufferInfo.height = extent.height;
            framebufferInfo.layers = 1;

            if (vkCreateFramebuffer(m_device->GetNativeDevice(), &framebufferInfo, nullptr, &m_framebuffers[i]) != VK_SUCCESS) {
                throw std::runtime_error("Failed to create framebuffer!");
            }
        }

        std::cout << "[VulkanFramebuffer] Created " << m_framebuffers.size() << " Framebuffers." << std::endl;
    }

    VulkanFramebuffer::~VulkanFramebuffer() {
        for (auto framebuffer : m_framebuffers) {
            vkDestroyFramebuffer(m_device->GetNativeDevice(), framebuffer, nullptr);
        }
    }

} // namespace Lizeral