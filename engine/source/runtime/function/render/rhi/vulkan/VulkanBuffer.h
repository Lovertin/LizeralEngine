#pragma once
#include <vulkan/vulkan.h>
#include <vma/vk_mem_alloc.h>

namespace Lizeral {

    class VulkanDevice;

    class VulkanBuffer {
    public:

        VulkanBuffer(VulkanDevice* device, 
                     VkDeviceSize size, 
                     VkBufferUsageFlags usage, 
                     VmaMemoryUsage memoryUsage);
        ~VulkanBuffer();

        // RAII
        VulkanBuffer(const VulkanBuffer&) = delete;
        VulkanBuffer& operator=(const VulkanBuffer&) = delete;

        VkBuffer GetNativeBuffer() const { return m_buffer; }
        VkDeviceSize GetSize() const { return m_size; }
        
        uint64_t GetDeviceAddress() const;

        void* Map();
        void Unmap();
        void WriteData(const void* data, VkDeviceSize size, VkDeviceSize offset = 0);

    private:
        VulkanDevice* m_device { nullptr };
        VkBuffer m_buffer { VK_NULL_HANDLE };
        VmaAllocation m_allocation { VK_NULL_HANDLE };
        VmaAllocationInfo m_allocInfo {};
        
        VkDeviceSize m_size { 0 };
        void* m_mappedData { nullptr };
    };

} // namespace Lizeral