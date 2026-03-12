#include "VulkanTLAS.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include <stdexcept>
#include <cstring>

namespace Lizeral {

    VulkanTLAS::VulkanTLAS(VulkanDevice* device, uint32_t maxFrames) : m_device(device), m_maxFrames(maxFrames) {
        // 初始化数组大小
        m_tlasBuffer.resize(m_maxFrames, VK_NULL_HANDLE);
        m_tlasMemory.resize(m_maxFrames, VK_NULL_HANDLE);
        m_tlasHandle.resize(m_maxFrames, VK_NULL_HANDLE);

        m_instanceBuffer.resize(m_maxFrames, VK_NULL_HANDLE);
        m_instanceMemory.resize(m_maxFrames, VK_NULL_HANDLE);
        m_instanceAddress.resize(m_maxFrames, 0);

        m_scratchBuffer.resize(m_maxFrames, VK_NULL_HANDLE);
        m_scratchMemory.resize(m_maxFrames, VK_NULL_HANDLE);

        m_maxInstanceCount.resize(m_maxFrames, 0);
        m_currentTlasSize.resize(m_maxFrames, 0);

        loadRTFunctions();
    }

    VulkanTLAS::~VulkanTLAS() {
        cleanup();
    }

    void VulkanTLAS::cleanup() {
        VkDevice logicalDevice = m_device->GetNativeDevice();

        // 循环清理每一帧的资源
        for (uint32_t i = 0; i < m_maxFrames; i++) {
            if (m_tlasHandle[i] != VK_NULL_HANDLE) pfn_vkDestroyAccelerationStructureKHR(logicalDevice, m_tlasHandle[i], nullptr);
            if (m_tlasBuffer[i] != VK_NULL_HANDLE) { vkDestroyBuffer(logicalDevice, m_tlasBuffer[i], nullptr); vkFreeMemory(logicalDevice, m_tlasMemory[i], nullptr); }
            if (m_instanceBuffer[i] != VK_NULL_HANDLE) { vkDestroyBuffer(logicalDevice, m_instanceBuffer[i], nullptr); vkFreeMemory(logicalDevice, m_instanceMemory[i], nullptr); }
            if (m_scratchBuffer[i] != VK_NULL_HANDLE) { vkDestroyBuffer(logicalDevice, m_scratchBuffer[i], nullptr); vkFreeMemory(logicalDevice, m_scratchMemory[i], nullptr); }
        }
    }

    void VulkanTLAS::Build(VkCommandBuffer cmd, uint32_t frameIndex, const std::vector<VkAccelerationStructureInstanceKHR>& instances) {
        if (instances.empty()) return;

        VkDevice logicalDevice = m_device->GetNativeDevice();
        uint32_t instanceCount = static_cast<uint32_t>(instances.size());
        
        // ★ 核心：计算当前应该操作哪一套缓冲 (0 或 1)
        uint32_t idx = frameIndex % m_maxFrames;

        // 1. 实例缓冲扩容逻辑
        if (instanceCount > m_maxInstanceCount[idx]) {
            if (m_instanceBuffer[idx] != VK_NULL_HANDLE) {
                vkDestroyBuffer(logicalDevice, m_instanceBuffer[idx], nullptr);
                vkFreeMemory(logicalDevice, m_instanceMemory[idx], nullptr);
            }
            m_maxInstanceCount[idx] = instanceCount * 2; // 预留一倍空间

            VkDeviceSize instanceBufferSize = sizeof(VkAccelerationStructureInstanceKHR) * m_maxInstanceCount[idx];
            allocateBuffer(
                instanceBufferSize, 
                VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT | VK_BUFFER_USAGE_ACCELERATION_STRUCTURE_BUILD_INPUT_READ_ONLY_BIT_KHR,
                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                m_instanceBuffer[idx], m_instanceMemory[idx], &m_instanceAddress[idx]
            );
        }

        // 2. 拷贝实例数据 (绝对安全，因为这块内存当前帧只有 CPU 在碰)
        void* mappedData;
        vkMapMemory(logicalDevice, m_instanceMemory[idx], 0, sizeof(VkAccelerationStructureInstanceKHR) * instanceCount, 0, &mappedData);
        memcpy(mappedData, instances.data(), sizeof(VkAccelerationStructureInstanceKHR) * instanceCount);
        vkUnmapMemory(logicalDevice, m_instanceMemory[idx]);

        // 3. 配置 Geometry 信息
        VkAccelerationStructureGeometryKHR geometry{};
        geometry.sType = VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_GEOMETRY_KHR;
        geometry.geometryType = VK_GEOMETRY_TYPE_INSTANCES_KHR;
        geometry.flags = VK_GEOMETRY_OPAQUE_BIT_KHR;
        geometry.geometry.instances.sType = VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_GEOMETRY_INSTANCES_DATA_KHR;
        geometry.geometry.instances.arrayOfPointers = VK_FALSE;
        geometry.geometry.instances.data.deviceAddress = m_instanceAddress[idx]; // ★ 用当前帧的地址

        // 4. 获取构建大小
        VkAccelerationStructureBuildGeometryInfoKHR buildInfo{};
        buildInfo.sType = VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_BUILD_GEOMETRY_INFO_KHR;
        buildInfo.type = VK_ACCELERATION_STRUCTURE_TYPE_TOP_LEVEL_KHR;
        buildInfo.flags = VK_BUILD_ACCELERATION_STRUCTURE_PREFER_FAST_TRACE_BIT_KHR;
        buildInfo.mode = VK_BUILD_ACCELERATION_STRUCTURE_MODE_BUILD_KHR;
        buildInfo.geometryCount = 1;
        buildInfo.pGeometries = &geometry;

        VkAccelerationStructureBuildSizesInfoKHR sizeInfo{VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_BUILD_SIZES_INFO_KHR};
        pfn_vkGetAccelerationStructureBuildSizesKHR(logicalDevice, VK_ACCELERATION_STRUCTURE_BUILD_TYPE_DEVICE_KHR, &buildInfo, &instanceCount, &sizeInfo);

        // 5. TLAS 本体扩容逻辑 (万一加了新物体，需要的树节点变多了)
        if (sizeInfo.accelerationStructureSize > m_currentTlasSize[idx]) {
            if (m_tlasHandle[idx] != VK_NULL_HANDLE) pfn_vkDestroyAccelerationStructureKHR(logicalDevice, m_tlasHandle[idx], nullptr);
            if (m_tlasBuffer[idx] != VK_NULL_HANDLE) { vkDestroyBuffer(logicalDevice, m_tlasBuffer[idx], nullptr); vkFreeMemory(logicalDevice, m_tlasMemory[idx], nullptr); }
            
            m_currentTlasSize[idx] = sizeInfo.accelerationStructureSize;

            allocateBuffer(
                sizeInfo.accelerationStructureSize,
                VK_BUFFER_USAGE_ACCELERATION_STRUCTURE_STORAGE_BIT_KHR | VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT,
                VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                m_tlasBuffer[idx], m_tlasMemory[idx]
            );

            VkAccelerationStructureCreateInfoKHR createInfo{VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_CREATE_INFO_KHR};
            createInfo.buffer = m_tlasBuffer[idx];
            createInfo.size = sizeInfo.accelerationStructureSize;
            createInfo.type = VK_ACCELERATION_STRUCTURE_TYPE_TOP_LEVEL_KHR;
            pfn_vkCreateAccelerationStructureKHR(logicalDevice, &createInfo, nullptr, &m_tlasHandle[idx]);
        }

        // 6. 草稿纸重分配 (安全，因为只摧毁当前帧旧的草稿纸，而上一帧的草稿纸在另一个 index 里完好无损)
        if (m_scratchBuffer[idx] != VK_NULL_HANDLE) {
            vkDestroyBuffer(logicalDevice, m_scratchBuffer[idx], nullptr);
            vkFreeMemory(logicalDevice, m_scratchMemory[idx], nullptr);
        }
        VkDeviceAddress scratchAddress;
        allocateBuffer(
            sizeInfo.buildScratchSize,
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
            m_scratchBuffer[idx], m_scratchMemory[idx], &scratchAddress
        );

        // 7. 录制命令
        buildInfo.dstAccelerationStructure = m_tlasHandle[idx]; // ★ 构建到当前帧的 TLAS 里
        buildInfo.scratchData.deviceAddress = scratchAddress;

        VkAccelerationStructureBuildRangeInfoKHR buildRange{};
        buildRange.primitiveCount = instanceCount;
        buildRange.primitiveOffset = 0;
        buildRange.firstVertex = 0;
        buildRange.transformOffset = 0;
        const VkAccelerationStructureBuildRangeInfoKHR* pBuildRange = &buildRange;

        pfn_vkCmdBuildAccelerationStructuresKHR(cmd, 1, &buildInfo, &pBuildRange);
    }

    // allocateBuffer 和 loadRTFunctions 保持原样不变 ...
    // (省略重复代码)
    void VulkanTLAS::allocateBuffer(VkDeviceSize size, VkBufferUsageFlags usage, VkMemoryPropertyFlags properties, VkBuffer& outBuffer, VkDeviceMemory& outMemory, VkDeviceAddress* outAddress) {
        VkDevice logicalDevice = m_device->GetNativeDevice();

        VkBufferCreateInfo bufferInfo{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        bufferInfo.size = size;
        bufferInfo.usage = usage;
        vkCreateBuffer(logicalDevice, &bufferInfo, nullptr, &outBuffer);

        VkMemoryRequirements memReqs;
        vkGetBufferMemoryRequirements(logicalDevice, outBuffer, &memReqs);

        VkPhysicalDeviceMemoryProperties memProps;
        vkGetPhysicalDeviceMemoryProperties(m_device->GetContext()->GetPhysicalDevice(), &memProps);
        uint32_t memTypeIndex = 0;
        for (uint32_t i = 0; i < memProps.memoryTypeCount; i++) {
            if ((memReqs.memoryTypeBits & (1 << i)) && (memProps.memoryTypes[i].propertyFlags & properties) == properties) {
                memTypeIndex = i; break;
            }
        }

        VkMemoryAllocateFlagsInfo allocFlagsInfo{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO};
        allocFlagsInfo.flags = VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT;
        VkMemoryAllocateInfo allocInfo{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        allocInfo.allocationSize = memReqs.size;
        allocInfo.memoryTypeIndex = memTypeIndex;
        allocInfo.pNext = &allocFlagsInfo;

        if (vkAllocateMemory(logicalDevice, &allocInfo, nullptr, &outMemory) != VK_SUCCESS) {
            throw std::runtime_error("Failed to allocate memory for TLAS!");
        }
        vkBindBufferMemory(logicalDevice, outBuffer, outMemory, 0);

        if (outAddress != nullptr) {
            VkBufferDeviceAddressInfo addressInfo{VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO};
            addressInfo.buffer = outBuffer;
            *outAddress = vkGetBufferDeviceAddress(logicalDevice, &addressInfo);
        }
    }

    void VulkanTLAS::loadRTFunctions() {
        VkDevice device = m_device->GetNativeDevice();
        pfn_vkGetAccelerationStructureBuildSizesKHR = (PFN_vkGetAccelerationStructureBuildSizesKHR)vkGetDeviceProcAddr(device, "vkGetAccelerationStructureBuildSizesKHR");
        pfn_vkCreateAccelerationStructureKHR = (PFN_vkCreateAccelerationStructureKHR)vkGetDeviceProcAddr(device, "vkCreateAccelerationStructureKHR");
        pfn_vkCmdBuildAccelerationStructuresKHR = (PFN_vkCmdBuildAccelerationStructuresKHR)vkGetDeviceProcAddr(device, "vkCmdBuildAccelerationStructuresKHR");
        pfn_vkDestroyAccelerationStructureKHR = (PFN_vkDestroyAccelerationStructureKHR)vkGetDeviceProcAddr(device, "vkDestroyAccelerationStructureKHR");
        pfn_vkGetAccelerationStructureDeviceAddressKHR = (PFN_vkGetAccelerationStructureDeviceAddressKHR)vkGetDeviceProcAddr(device, "vkGetAccelerationStructureDeviceAddressKHR");
    }

}