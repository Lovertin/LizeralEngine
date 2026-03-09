#include "VulkanCommandBuffer.h"
#include "VulkanDevice.h"
#include "VulkanCommandPool.h"
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanCommandBuffer::VulkanCommandBuffer(VulkanDevice* device, VulkanCommandPool* pool) 
        : m_device(device), m_pool(pool) {
        
        VkCommandBufferAllocateInfo allocInfo{};
        allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
        allocInfo.commandPool = m_pool->GetNativePool();
        allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY; // 可以直接提交给队列的主级缓冲
        allocInfo.commandBufferCount = 1;

        if (vkAllocateCommandBuffers(m_device->GetNativeDevice(), &allocInfo, &m_commandBuffer) != VK_SUCCESS) {
            throw std::runtime_error("Failed to allocate command buffer!");
        }
    }

    VulkanCommandBuffer::~VulkanCommandBuffer() {
        if (m_commandBuffer != VK_NULL_HANDLE) {
            // 归还清单给命令池
            vkFreeCommandBuffers(m_device->GetNativeDevice(), m_pool->GetNativePool(), 1, &m_commandBuffer);
        }
    }

    void VulkanCommandBuffer::Begin(VkCommandBufferUsageFlags flags) {
        VkCommandBufferBeginInfo beginInfo{};
        beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        beginInfo.flags = flags; // 比如传入 VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT

        if (vkBeginCommandBuffer(m_commandBuffer, &beginInfo) != VK_SUCCESS) {
            throw std::runtime_error("Failed to begin recording command buffer!");
        }
    }

    void VulkanCommandBuffer::End() {
        if (vkEndCommandBuffer(m_commandBuffer) != VK_SUCCESS) {
            throw std::runtime_error("Failed to record command buffer!");
        }
    }

    void VulkanCommandBuffer::SubmitAndIdle() {
        VkSubmitInfo submitInfo{};
        submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &m_commandBuffer;

        // 1. 把清单扔给传送带 (Graphics Queue)
        vkQueueSubmit(m_device->GetGraphicsQueue(), 1, &submitInfo, VK_NULL_HANDLE);

        // 2. 暴力阻塞 CPU，强行等待这条传送带停下来 (专门用于初始化阶段的同步)
        vkQueueWaitIdle(m_device->GetGraphicsQueue());
    }

} // namespace Lizeral