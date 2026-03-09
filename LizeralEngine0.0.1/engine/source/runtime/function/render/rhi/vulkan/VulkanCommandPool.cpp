#include "VulkanCommandPool.h"
#include "VulkanDevice.h"
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanCommandPool::VulkanCommandPool(VulkanDevice* device) : m_device(device) {
        VkCommandPoolCreateInfo poolInfo{};
        poolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
        
        // 【关键标志】：允许我们单独重置池子里的某一张命令缓冲（而不是每次都要销毁整个池子）
        poolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT; 
        
        // 绑定到我们在 Device 中找到的图形队列家族 (Index 0)
        poolInfo.queueFamilyIndex = m_device->GetGraphicsQueueFamily();

        if (vkCreateCommandPool(m_device->GetNativeDevice(), &poolInfo, nullptr, &m_commandPool) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create command pool!");
        }
        
        std::cout << "[VulkanCommandPool] Command Pool created successfully." << std::endl;
    }

    VulkanCommandPool::~VulkanCommandPool() {
        if (m_commandPool != VK_NULL_HANDLE) {
            vkDestroyCommandPool(m_device->GetNativeDevice(), m_commandPool, nullptr);
            m_commandPool = VK_NULL_HANDLE;
        }
    }

} // namespace Lizeral