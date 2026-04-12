#pragma once
#include <cstdint>
#include <memory>
#include <vector>
#include "runtime/function/render/rhi/vulkan/VulkanBuffer.h"
#include "runtime/resource/material/materialAsset.h"
#include "runtime/resource/modelResourceData.h"
#include "runtime/resource/meshletAssetTypes.h"

namespace Lizeral {

    class VulkanBLAS;

    struct VulkanModelResource {
        std::unique_ptr<VulkanBuffer> vertexBuffer;
        std::unique_ptr<VulkanBuffer> meshletBuffer;
        std::unique_ptr<VulkanBuffer> indexBuffer;
        std::unique_ptr<VulkanBuffer> boundsBuffer;
        std::unique_ptr<VulkanBuffer> materialBuffer;
        std::unique_ptr<VulkanBuffer> globalIndexBuffer;
        std::unique_ptr<VulkanBuffer> primitiveMaterialIdBuffer;

        uint64_t vAddr = 0;         
        uint64_t mAddr = 0;         
        uint64_t iAddr = 0;         
        uint64_t bAddr = 0;         
        uint64_t matAddr = 0;       
        uint64_t globalIAddr = 0;
        uint64_t primMatIdAddr = 0;
        
        uint32_t totalMeshlets = 0; 
        uint32_t textureOffset = 0; 
        uint32_t textureCount = 0;  
        uint32_t materialCount = 0;
        std::vector<MaterialData> materialsCpu;
        std::vector<Resource::MaterialAsset> materialAssets;
        std::vector<Resource::RuntimeMeshAssetData> meshAssets;
        std::shared_ptr<VulkanBLAS> blas;

        bool IsValid() const { return vAddr != 0; }
    };
}
