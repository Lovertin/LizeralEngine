#pragma once
#include <vulkan/vulkan.h>

namespace Lizeral {

    class VulkanDevice; 

    class VulkanCommandPool {
    public:
        VulkanCommandPool(VulkanDevice* device);
        ~VulkanCommandPool();

        VulkanCommandPool(const VulkanCommandPool&) = delete;
        VulkanCommandPool& operator=(const VulkanCommandPool&) = delete;

        VkCommandPool GetNativePool() const { return m_commandPool; }

    private:
        VulkanDevice* m_device { nullptr };
        VkCommandPool m_commandPool { VK_NULL_HANDLE };
    };

} // namespace Lizeral