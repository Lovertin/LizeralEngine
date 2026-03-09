#pragma once
#include <vulkan/vulkan.h>

namespace Lizeral {

    class VulkanDevice; // 前向声明

    class VulkanCommandPool {
    public:
        // 构造函数需要知道给哪个设备、哪个车间（队列）建池子
        VulkanCommandPool(VulkanDevice* device);
        ~VulkanCommandPool();

        // 禁用拷贝
        VulkanCommandPool(const VulkanCommandPool&) = delete;
        VulkanCommandPool& operator=(const VulkanCommandPool&) = delete;

        VkCommandPool GetNativePool() const { return m_commandPool; }

    private:
        VulkanDevice* m_device { nullptr };
        VkCommandPool m_commandPool { VK_NULL_HANDLE };
    };

} // namespace Lizeral