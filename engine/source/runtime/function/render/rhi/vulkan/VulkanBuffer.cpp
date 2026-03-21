#include "VulkanBuffer.h"
#include "VulkanDevice.h"
#include <stdexcept>
#include <cstring>

namespace Lizeral {

    VulkanBuffer::VulkanBuffer(VulkanDevice* device, VkDeviceSize size, VkBufferUsageFlags usage, VmaMemoryUsage memoryUsage)
        : m_device(device), m_size(size) {
        
        VkBufferCreateInfo bufferInfo{};
        bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
        bufferInfo.size = size;
        bufferInfo.usage = usage;
        bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

        VmaAllocationCreateInfo allocInfo{};
        allocInfo.usage = memoryUsage;
        if (memoryUsage == VMA_MEMORY_USAGE_CPU_TO_GPU || memoryUsage == VMA_MEMORY_USAGE_CPU_ONLY) {
            allocInfo.flags = VMA_ALLOCATION_CREATE_MAPPED_BIT; 
        }

        if (vmaCreateBuffer(m_device->GetAllocator(), &bufferInfo, &allocInfo, 
                            &m_buffer, &m_allocation, &m_allocInfo) != VK_SUCCESS) {
            throw std::runtime_error("Failed to allocate GPU Buffer via VMA!");
        }

        m_mappedData = m_allocInfo.pMappedData; 
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
        
        return vkGetBufferDeviceAddress(m_device->GetNativeDevice(), &addressInfo);
    }

    void* VulkanBuffer::Map() {
        if (m_allocInfo.pMappedData) return m_allocInfo.pMappedData;
    
        if (!m_mappedData) vmaMapMemory(m_device->GetAllocator(), m_allocation, &m_mappedData);
        return m_mappedData;
    }

    void VulkanBuffer::Unmap() {
        if (m_allocInfo.pMappedData) return; 
    
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
        }
    }

} // namespace Lizeral