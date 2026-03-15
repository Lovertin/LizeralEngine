#pragma once
#include <cstdint>
#include <memory>
#include "runtime/function/render/rhi/vulkan/VulkanBuffer.h"

namespace Lizeral {

    class VulkanBLAS;

    struct VulkanModelResource {
        // ★ 让 Resource 直接拥有这些 Buffer 的生命周期
        std::unique_ptr<VulkanBuffer> vertexBuffer;
        std::unique_ptr<VulkanBuffer> meshletBuffer;
        std::unique_ptr<VulkanBuffer> indexBuffer;
        std::unique_ptr<VulkanBuffer> boundsBuffer;
        std::unique_ptr<VulkanBuffer> materialBuffer;
        std::unique_ptr<VulkanBuffer> globalIndexBuffer;

        // 缓存 BDA 地址，方便每帧极速组装 PushConstants
        uint64_t vAddr = 0;         
        uint64_t mAddr = 0;         
        uint64_t iAddr = 0;         
        uint64_t bAddr = 0;         
        uint64_t matAddr = 0;       
        uint64_t globalIAddr = 0;
        
        uint32_t totalMeshlets = 0; 
        uint32_t textureOffset = 0; 
        uint32_t textureCount = 0;  
        std::shared_ptr<VulkanBLAS> blas;

        bool IsValid() const { return vAddr != 0; }
    };
}