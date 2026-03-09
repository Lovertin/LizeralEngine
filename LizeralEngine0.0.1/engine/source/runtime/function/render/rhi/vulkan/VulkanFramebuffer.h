#pragma once
#include <vulkan/vulkan.h>
#include <vector>

namespace Lizeral {

    class VulkanDevice;
    class VulkanRenderPass;
    class VulkanSwapchain;

    class VulkanFramebuffer {
    public:
        // 构造函数：把 Device、说明书(RenderPass) 和 相框集合(Swapchain) 传进来
        VulkanFramebuffer(VulkanDevice* device, VulkanRenderPass* renderPass, VulkanSwapchain* swapchain);
        ~VulkanFramebuffer();

        // 禁用拷贝
        VulkanFramebuffer(const VulkanFramebuffer&) = delete;
        VulkanFramebuffer& operator=(const VulkanFramebuffer&) = delete;

        // 获取所有的实体画板
        const std::vector<VkFramebuffer>& GetNativeFramebuffers() const { return m_framebuffers; }

    private:
        VulkanDevice* m_device { nullptr };
        std::vector<VkFramebuffer> m_framebuffers;
    };

} // namespace Lizeral