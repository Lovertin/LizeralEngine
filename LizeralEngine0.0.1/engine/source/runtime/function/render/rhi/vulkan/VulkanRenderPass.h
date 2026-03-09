#pragma once
#include <vulkan/vulkan.h>

namespace Lizeral {

    class VulkanDevice;
    class VulkanSwapchain;

    class VulkanRenderPass {
    public:
        // 构造时需要 Device(负责创建) 和 Swapchain(为了知道图片是什么颜色格式)
        VulkanRenderPass(VulkanDevice* device, VulkanSwapchain* swapchain);
        ~VulkanRenderPass();

        VulkanRenderPass(const VulkanRenderPass&) = delete;
        VulkanRenderPass& operator=(const VulkanRenderPass&) = delete;

        VkRenderPass GetNativeRenderPass() const { return m_renderPass; }

    private:
        VulkanDevice* m_device { nullptr };
        VkRenderPass m_renderPass { VK_NULL_HANDLE };
    };

} // namespace Lizeral