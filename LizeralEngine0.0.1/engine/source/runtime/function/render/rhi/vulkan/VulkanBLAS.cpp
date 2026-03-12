#include "VulkanBLAS.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandPool.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandBuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanBLAS::VulkanBLAS(VulkanDevice* device, VulkanCommandPool* commandPool, 
                           VkDeviceAddress vertexAddress, uint32_t vertexCount, uint32_t vertexStride,
                           VkDeviceAddress indexAddress,  uint32_t indexCount)
        : m_device(device) {
        
        loadRTFunctions();

        VkDevice logicalDevice = m_device->GetNativeDevice();

        // 1. 描述几何体 (告诉 GPU 我们的模型长什么样)
        VkAccelerationStructureGeometryKHR geometry{};
        geometry.sType = VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_GEOMETRY_KHR;
        geometry.flags = VK_GEOMETRY_OPAQUE_BIT_KHR; // 不透明，加速光线求交
        geometry.geometryType = VK_GEOMETRY_TYPE_TRIANGLES_KHR;

        geometry.geometry.triangles.sType = VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_GEOMETRY_TRIANGLES_DATA_KHR;
        geometry.geometry.triangles.vertexFormat = VK_FORMAT_R32G32B32_SFLOAT; // 假设顶点位置是 vec3
        geometry.geometry.triangles.vertexData.deviceAddress = vertexAddress;
        geometry.geometry.triangles.maxVertex = vertexCount - 1;
        geometry.geometry.triangles.vertexStride = vertexStride;
        // 你用的是 MicroIndices (可能是 uint8 或 uint32，假设你打包成了 uint32，请根据 Meshlet 实际格式修改)
        geometry.geometry.triangles.indexType = VK_INDEX_TYPE_UINT32; 
        geometry.geometry.triangles.indexData.deviceAddress = indexAddress;

        // 2. 询问 GPU 需要多大的显存来装这个 BLAS 和 临时计算空间 (Scratch)
        VkAccelerationStructureBuildGeometryInfoKHR buildInfo{};
        buildInfo.sType = VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_BUILD_GEOMETRY_INFO_KHR;
        buildInfo.type = VK_ACCELERATION_STRUCTURE_TYPE_BOTTOM_LEVEL_KHR;
        buildInfo.flags = VK_BUILD_ACCELERATION_STRUCTURE_PREFER_FAST_TRACE_BIT_KHR;
        buildInfo.mode = VK_BUILD_ACCELERATION_STRUCTURE_MODE_BUILD_KHR;
        buildInfo.geometryCount = 1;
        buildInfo.pGeometries = &geometry;

        uint32_t numTriangles = indexCount / 3;

        VkAccelerationStructureBuildSizesInfoKHR sizeInfo{VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_BUILD_SIZES_INFO_KHR};
        pfn_vkGetAccelerationStructureBuildSizesKHR(logicalDevice, VK_ACCELERATION_STRUCTURE_BUILD_TYPE_DEVICE_KHR, &buildInfo, &numTriangles, &sizeInfo);

        // 3. 分配 BLAS 专用的显存 Buffer
        VkBufferCreateInfo bufferInfo{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        bufferInfo.size = sizeInfo.accelerationStructureSize;
        bufferInfo.usage = VK_BUFFER_USAGE_ACCELERATION_STRUCTURE_STORAGE_BIT_KHR | VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
        vkCreateBuffer(logicalDevice, &bufferInfo, nullptr, &m_asBuffer);

        VkMemoryRequirements memReqs;
        vkGetBufferMemoryRequirements(logicalDevice, m_asBuffer, &memReqs);

        // (此处简化内存分配，你可以复用 sandbox 里的内存查找逻辑，或者直接用 VMA)
        VkPhysicalDeviceMemoryProperties memProps;
        vkGetPhysicalDeviceMemoryProperties(m_device->GetContext()->GetPhysicalDevice(), &memProps);
        uint32_t memTypeIndex = 0;
        for (uint32_t i = 0; i < memProps.memoryTypeCount; i++) {
            if ((memReqs.memoryTypeBits & (1 << i)) && 
                (memProps.memoryTypes[i].propertyFlags & VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT) == VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT) {
                memTypeIndex = i; break;
            }
        }

        VkMemoryAllocateFlagsInfo allocFlagsInfo{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO};
        allocFlagsInfo.flags = VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT;
        VkMemoryAllocateInfo allocInfo{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        allocInfo.allocationSize = memReqs.size;
        allocInfo.memoryTypeIndex = memTypeIndex;
        allocInfo.pNext = &allocFlagsInfo;
        vkAllocateMemory(logicalDevice, &allocInfo, nullptr, &m_asMemory);
        vkBindBufferMemory(logicalDevice, m_asBuffer, m_asMemory, 0);

        // 4. 正式创建 AS 对象
        VkAccelerationStructureCreateInfoKHR createInfo{VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_CREATE_INFO_KHR};
        createInfo.buffer = m_asBuffer;
        createInfo.size = sizeInfo.accelerationStructureSize;
        createInfo.type = VK_ACCELERATION_STRUCTURE_TYPE_BOTTOM_LEVEL_KHR;
        pfn_vkCreateAccelerationStructureKHR(logicalDevice, &createInfo, nullptr, &m_blasHandle);

        // 获取 AS 的设备地址 (TLAS 构建时需要用到)
        VkAccelerationStructureDeviceAddressInfoKHR addressInfo{VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_DEVICE_ADDRESS_INFO_KHR};
        addressInfo.accelerationStructure = m_blasHandle;
        m_blasAddress = pfn_vkGetAccelerationStructureDeviceAddressKHR(logicalDevice, &addressInfo);

        // 5. 分配一个临时 Scratch Buffer (GPU 构建 AS 时的草稿纸)
        VkBuffer scratchBuffer;
        VkDeviceMemory scratchMemory;
        bufferInfo.size = sizeInfo.buildScratchSize;
        bufferInfo.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
        vkCreateBuffer(logicalDevice, &bufferInfo, nullptr, &scratchBuffer);
        
        vkGetBufferMemoryRequirements(logicalDevice, scratchBuffer, &memReqs);
        allocInfo.allocationSize = memReqs.size;
        // 简单复用刚才的 memTypeIndex (要求 Device Local)
        vkAllocateMemory(logicalDevice, &allocInfo, nullptr, &scratchMemory);
        vkBindBufferMemory(logicalDevice, scratchBuffer, scratchMemory, 0);

        VkBufferDeviceAddressInfo scratchAddressInfo{VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO};
        scratchAddressInfo.buffer = scratchBuffer;
        VkDeviceAddress scratchAddress = vkGetBufferDeviceAddress(logicalDevice, &scratchAddressInfo);

        // 6. 录制命令，让 GPU 真正开始“干活”构建树结构
        VulkanCommandBuffer cmdBuf(m_device, commandPool);
        cmdBuf.Begin(VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT);
        VkCommandBuffer cmd = cmdBuf.GetNativeBuffer();

        buildInfo.dstAccelerationStructure = m_blasHandle;
        buildInfo.scratchData.deviceAddress = scratchAddress;

        VkAccelerationStructureBuildRangeInfoKHR buildRange{};
        buildRange.primitiveCount = numTriangles;
        buildRange.primitiveOffset = 0;
        buildRange.firstVertex = 0;
        buildRange.transformOffset = 0;
        const VkAccelerationStructureBuildRangeInfoKHR* pBuildRange = &buildRange;

        pfn_vkCmdBuildAccelerationStructuresKHR(cmd, 1, &buildInfo, &pBuildRange);

        cmdBuf.End();
        cmdBuf.SubmitAndIdle(); // 阻塞等待 GPU 构建完毕

        // 构建完了，草稿纸可以扔掉了
        vkDestroyBuffer(logicalDevice, scratchBuffer, nullptr);
        vkFreeMemory(logicalDevice, scratchMemory, nullptr);

        std::cout << "[VulkanBLAS] Successfully built BLAS for model. Triangles: " << numTriangles << std::endl;
    }

    VulkanBLAS::~VulkanBLAS() {
        if (m_blasHandle != VK_NULL_HANDLE) {
            pfn_vkDestroyAccelerationStructureKHR(m_device->GetNativeDevice(), m_blasHandle, nullptr);
        }
        if (m_asBuffer != VK_NULL_HANDLE) {
            vkDestroyBuffer(m_device->GetNativeDevice(), m_asBuffer, nullptr);
            vkFreeMemory(m_device->GetNativeDevice(), m_asMemory, nullptr);
        }
    }

    void VulkanBLAS::loadRTFunctions() {
        VkDevice device = m_device->GetNativeDevice();
        pfn_vkGetAccelerationStructureBuildSizesKHR = (PFN_vkGetAccelerationStructureBuildSizesKHR)vkGetDeviceProcAddr(device, "vkGetAccelerationStructureBuildSizesKHR");
        pfn_vkCreateAccelerationStructureKHR = (PFN_vkCreateAccelerationStructureKHR)vkGetDeviceProcAddr(device, "vkCreateAccelerationStructureKHR");
        pfn_vkCmdBuildAccelerationStructuresKHR = (PFN_vkCmdBuildAccelerationStructuresKHR)vkGetDeviceProcAddr(device, "vkCmdBuildAccelerationStructuresKHR");
        pfn_vkDestroyAccelerationStructureKHR = (PFN_vkDestroyAccelerationStructureKHR)vkGetDeviceProcAddr(device, "vkDestroyAccelerationStructureKHR");
        pfn_vkGetAccelerationStructureDeviceAddressKHR = (PFN_vkGetAccelerationStructureDeviceAddressKHR)vkGetDeviceProcAddr(device, "vkGetAccelerationStructureDeviceAddressKHR");
        
        if (!pfn_vkCreateAccelerationStructureKHR) {
            throw std::runtime_error("Failed to load Vulkan Ray Tracing function pointers!");
        }
    }

}