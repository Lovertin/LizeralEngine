#pragma once
#include <cstdint>

namespace Lizeral {
    // 纯数据结构，没有任何 Vulkan/OpenGL 的杂质
    struct VulkanModelResource {
        uint64_t vAddr = 0;         // 顶点缓冲 BDA
        uint64_t mAddr = 0;         // Meshlet 缓冲 BDA
        uint64_t iAddr = 0;         // 微索引缓冲 BDA
        uint64_t bAddr = 0;         // 包围盒缓冲 BDA
        uint64_t matAddr = 0;       // PBR 材质缓冲 BDA
        uint32_t totalMeshlets = 0; // Meshlet 总数
        
        // 记录一下这个模型占用了全局贴图池的多少张图，方便以后清理
        uint32_t textureOffset = 0; 
        uint32_t textureCount = 0;  

        bool IsValid() const { return vAddr != 0; }
    };
}