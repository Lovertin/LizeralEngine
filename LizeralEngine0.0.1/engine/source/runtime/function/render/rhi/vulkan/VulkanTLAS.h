#pragma once
#include <vulkan/vulkan.h>
#include <vector>

namespace Lizeral {

    class VulkanDevice;
    class VulkanCommandPool;

    class VulkanTLAS {
    public:
        VulkanTLAS(VulkanDevice* device);
        ~VulkanTLAS();

        // 禁用拷贝
        VulkanTLAS(const VulkanTLAS&) = delete;
        VulkanTLAS& operator=(const VulkanTLAS&) = delete;

        // 核心方法：传入一组实例，构建 TLAS
        // cmd: 用于录制构建命令的 CommandBuffer
        void Build(VkCommandBuffer cmd, const std::vector<VkAccelerationStructureInstanceKHR>& instances);

        VkAccelerationStructureKHR GetHandle() const { return m_tlasHandle; }

    private:
        VulkanDevice* m_device{ nullptr };

        // TLAS 显存缓冲
        VkBuffer m_tlasBuffer{ VK_NULL_HANDLE };
        VkDeviceMemory m_tlasMemory{ VK_NULL_HANDLE };
        VkAccelerationStructureKHR m_tlasHandle{ VK_NULL_HANDLE };

        // 实例数据的缓冲 (存在 GPU 上，供 TLAS 读取)
        VkBuffer m_instanceBuffer{ VK_NULL_HANDLE };
        VkDeviceMemory m_instanceMemory{ VK_NULL_HANDLE };
        VkDeviceAddress m_instanceAddress{ 0 };

        // 草稿纸缓冲
        VkBuffer m_scratchBuffer{ VK_NULL_HANDLE };
        VkDeviceMemory m_scratchMemory{ VK_NULL_HANDLE };

        uint32_t m_maxInstanceCount{ 0 }; // 记录当前缓冲容量，避免频繁重新分配

        void loadRTFunctions();
        void allocateBuffer(VkDeviceSize size, VkBufferUsageFlags usage, VkMemoryPropertyFlags properties, VkBuffer& outBuffer, VkDeviceMemory& outMemory, VkDeviceAddress* outAddress = nullptr);
        void cleanup();

        // 函数指针
        PFN_vkGetAccelerationStructureBuildSizesKHR pfn_vkGetAccelerationStructureBuildSizesKHR = nullptr;
        PFN_vkCreateAccelerationStructureKHR pfn_vkCreateAccelerationStructureKHR = nullptr;
        PFN_vkCmdBuildAccelerationStructuresKHR pfn_vkCmdBuildAccelerationStructuresKHR = nullptr;
        PFN_vkDestroyAccelerationStructureKHR pfn_vkDestroyAccelerationStructureKHR = nullptr;
        PFN_vkGetAccelerationStructureDeviceAddressKHR pfn_vkGetAccelerationStructureDeviceAddressKHR = nullptr;
    };

}