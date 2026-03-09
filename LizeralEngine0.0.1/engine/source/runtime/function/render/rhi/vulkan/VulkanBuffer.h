#pragma once
#include <vulkan/vulkan.h>
#include <vma/vk_mem_alloc.h>

namespace Lizeral {

    class VulkanDevice; // 前向声明

    class VulkanBuffer {
    public:
        // 构造函数：需要我们的逻辑设备、所需大小、用途(比如是顶点还是SSBO)、以及内存域(VMA标志)
        VulkanBuffer(VulkanDevice* device, 
                     VkDeviceSize size, 
                     VkBufferUsageFlags usage, 
                     VmaMemoryUsage memoryUsage);
        ~VulkanBuffer();

        // 禁用拷贝，防止同一个显存被多次释放 (RAII 原则)
        VulkanBuffer(const VulkanBuffer&) = delete;
        VulkanBuffer& operator=(const VulkanBuffer&) = delete;

        // --- 核心获取接口 ---
        VkBuffer GetNativeBuffer() const { return m_buffer; }
        VkDeviceSize GetSize() const { return m_size; }
        
        // 【关键！】获取 BDA (64位显存指针)，这是未来做 Mesh Shader 和 Nanite 的入场券！
        uint64_t GetDeviceAddress() const;

        // --- CPU 映射接口 (仅对 CPU_TO_GPU 内存有效) ---
        void* Map();
        void Unmap();
        // 快速写入辅助函数
        void WriteData(const void* data, VkDeviceSize size, VkDeviceSize offset = 0);

    private:
        VulkanDevice* m_device { nullptr };
        VkBuffer m_buffer { VK_NULL_HANDLE };
        VmaAllocation m_allocation { VK_NULL_HANDLE };
        VmaAllocationInfo m_allocInfo {}; // 记录分配详情
        
        VkDeviceSize m_size { 0 };
        void* m_mappedData { nullptr };
    };

} // namespace Lizeral