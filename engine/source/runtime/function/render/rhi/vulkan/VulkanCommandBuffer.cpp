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
        allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        allocInfo.commandBufferCount = 1;

        if (vkAllocateCommandBuffers(m_device->GetNativeDevice(), &allocInfo, &m_commandBuffer) != VK_SUCCESS) {
            throw std::runtime_error("Failed to allocate command buffer!");
        }
    }

    VulkanCommandBuffer::~VulkanCommandBuffer() {
        if (m_commandBuffer != VK_NULL_HANDLE) {

            vkFreeCommandBuffers(m_device->GetNativeDevice(), m_pool->GetNativePool(), 1, &m_commandBuffer);
        }
    }

    void VulkanCommandBuffer::Begin(VkCommandBufferUsageFlags flags) {
        VkCommandBufferBeginInfo beginInfo{};
        beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        beginInfo.flags = flags; // VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT

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

        vkQueueSubmit(m_device->GetGraphicsQueue(), 1, &submitInfo, VK_NULL_HANDLE);

        vkQueueWaitIdle(m_device->GetGraphicsQueue());
    }

} // namespace Lizeral