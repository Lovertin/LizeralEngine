#pragma once
#include <vulkan/vulkan.h>
#include <memory>
#include <vector>

namespace Lizeral {

    class VulkanDevice;
    class VulkanCommandPool;
    class VulkanBuffer;

    class VulkanBLAS {
    public:
        VulkanBLAS(VulkanDevice* device, VulkanCommandPool* commandPool, 
                   VkDeviceAddress vertexAddress, uint32_t vertexCount, uint32_t vertexStride,
                   VkDeviceAddress indexAddress,  uint32_t indexCount);
        ~VulkanBLAS();

        VulkanBLAS(const VulkanBLAS&) = delete;
        VulkanBLAS& operator=(const VulkanBLAS&) = delete;

        VkAccelerationStructureKHR GetHandle() const { return m_blasHandle; }
        VkDeviceAddress GetDeviceAddress() const { return m_blasAddress; }

    private:
        VulkanDevice* m_device{ nullptr };
        
        VkBuffer m_asBuffer{ VK_NULL_HANDLE };
        VkDeviceMemory m_asMemory{ VK_NULL_HANDLE }; 
        
        VkAccelerationStructureKHR m_blasHandle{ VK_NULL_HANDLE };
        VkDeviceAddress m_blasAddress{ 0 };

        void loadRTFunctions();

        PFN_vkGetAccelerationStructureBuildSizesKHR pfn_vkGetAccelerationStructureBuildSizesKHR = nullptr;
        PFN_vkCreateAccelerationStructureKHR pfn_vkCreateAccelerationStructureKHR = nullptr;
        PFN_vkCmdBuildAccelerationStructuresKHR pfn_vkCmdBuildAccelerationStructuresKHR = nullptr;
        PFN_vkDestroyAccelerationStructureKHR pfn_vkDestroyAccelerationStructureKHR = nullptr;
        PFN_vkGetAccelerationStructureDeviceAddressKHR pfn_vkGetAccelerationStructureDeviceAddressKHR = nullptr;
    };

}