#pragma once
#include <vulkan/vulkan.h>
#include <vector>

namespace Lizeral {

    class VulkanDevice;

    class VulkanTLAS {
    public:
        // 增加 maxFrames 参数，默认 2（对应 Ping-Pong）
        VulkanTLAS(VulkanDevice* device, uint32_t maxFrames = 2);
        ~VulkanTLAS();

        // 禁用拷贝
        VulkanTLAS(const VulkanTLAS&) = delete;
        VulkanTLAS& operator=(const VulkanTLAS&) = delete;

        // ★ 核心改变：传入当前的 frameIndex，决定使用哪一套缓冲
        void Build(VkCommandBuffer cmd, uint32_t frameIndex, const std::vector<VkAccelerationStructureInstanceKHR>& instances);

        VkAccelerationStructureKHR GetHandle(uint32_t frameIndex) const { 
            return m_tlasHandle[frameIndex % m_maxFrames]; 
        }

    private:
        VulkanDevice* m_device{ nullptr };
        uint32_t m_maxFrames{ 2 };

        // ★ 全部变成数组，支持多重缓冲
        std::vector<VkBuffer> m_tlasBuffer;
        std::vector<VkDeviceMemory> m_tlasMemory;
        std::vector<VkAccelerationStructureKHR> m_tlasHandle;

        std::vector<VkBuffer> m_instanceBuffer;
        std::vector<VkDeviceMemory> m_instanceMemory;
        std::vector<VkDeviceAddress> m_instanceAddress;

        std::vector<VkBuffer> m_scratchBuffer;
        std::vector<VkDeviceMemory> m_scratchMemory;

        // 记录每个帧缓冲的当前容量
        std::vector<uint32_t> m_maxInstanceCount; 
        std::vector<VkDeviceSize> m_currentTlasSize; 

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