#include "VulkanBuffer.h"
#include "VulkanDevice.h"
#include <stdexcept>
#include <cstring>

namespace Lizeral {

    VulkanBuffer::VulkanBuffer(VulkanDevice* device, VkDeviceSize size, VkBufferUsageFlags usage, VmaMemoryUsage memoryUsage)
        : m_device(device), m_size(size) {
        
        // 1. 定义 Buffer 的基本信息
        VkBufferCreateInfo bufferInfo{};
        bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
        bufferInfo.size = size;
        bufferInfo.usage = usage;
        bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE; // 默认仅被一个队列家族拥有

        // 2. 定义 VMA 怎么去系统里拿内存
        VmaAllocationCreateInfo allocInfo{};
        allocInfo.usage = memoryUsage;
        // 如果是经常需要 CPU 映射的，加上随机访问标志
        if (memoryUsage == VMA_MEMORY_USAGE_CPU_TO_GPU || memoryUsage == VMA_MEMORY_USAGE_CPU_ONLY) {
            allocInfo.flags = VMA_ALLOCATION_CREATE_MAPPED_BIT; // VMA 优化：创建时直接映射好
        }

        // 3. 召唤 VMA 一键分配 Buffer + 绑定内存！
        if (vmaCreateBuffer(m_device->GetAllocator(), &bufferInfo, &allocInfo, 
                            &m_buffer, &m_allocation, &m_allocInfo) != VK_SUCCESS) {
            throw std::runtime_error("Failed to allocate GPU Buffer via VMA!");
        }

        m_mappedData = m_allocInfo.pMappedData; // 如果上面加了 MAPPED_BIT，这里直接就能拿到指针
    }

    VulkanBuffer::~VulkanBuffer() {
        if (m_buffer != VK_NULL_HANDLE && m_allocation != VK_NULL_HANDLE) {
            vmaDestroyBuffer(m_device->GetAllocator(), m_buffer, m_allocation);
        }
    }

    uint64_t VulkanBuffer::GetDeviceAddress() const {
        VkBufferDeviceAddressInfo addressInfo{};
        addressInfo.sType = VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO;
        addressInfo.buffer = m_buffer;
        
        // 注意：要在未来启用这个功能，我们需要在创建逻辑设备时开启 BDA 特性
        return vkGetBufferDeviceAddress(m_device->GetNativeDevice(), &addressInfo);
    }

    void* VulkanBuffer::Map() {
        if (!m_mappedData) {
            vmaMapMemory(m_device->GetAllocator(), m_allocation, &m_mappedData);
        }
        return m_mappedData;
    }

    void VulkanBuffer::Unmap() {
        if (m_mappedData) {
            vmaUnmapMemory(m_device->GetAllocator(), m_allocation);
            m_mappedData = nullptr;
        }
    }

    void VulkanBuffer::WriteData(const void* data, VkDeviceSize size, VkDeviceSize offset) {
        if (size > m_size) throw std::runtime_error("Buffer overflow!");
        
        void* mapped = Map();
        if (mapped) {
            memcpy(static_cast<char*>(mapped) + offset, data, size);
            // 如果内存不是 HOST_COHERENT 的，这里理论上需要调用 vmaFlushAllocation
            // 但为了简化，我们假设常用的 Uniform Buffer 都配置了 Coherent。
        }
    }

} // namespace Lizeral