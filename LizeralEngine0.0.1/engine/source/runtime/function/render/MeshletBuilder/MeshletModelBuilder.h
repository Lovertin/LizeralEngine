#pragma once

#include <vector>
#include <string>
#include <cstdint>

namespace Lizeral {

    // 1. 标准化的顶点结构 (内存必须严格对齐，以便 GPU 抓取)
    struct MeshletVertex {
        float pos[3];
        float normal[3];
        float uv[2];
    };

    // 2. Meshlet 描述符 (每个实例代表一个最多 64 顶点的微集群)
    struct MeshletDescriptor {
        uint32_t vertexOffset;    // 这个集群的顶点在【全局顶点数组】里的起始偏移
        uint32_t triangleOffset;  // 这个集群的索引在【全局微型索引数组】里的起始偏移
        uint32_t vertexCount;     // 顶点数量 (最多 64)
        uint32_t triangleCount;   // 三角形数量 (最多 126)
    };

    // 3. 构建器类：将 Assimp 和 MeshOptimizer 结合的终极兵工厂
    class MeshletModelBuilder {
    public:
        // 核心执行函数：加载模型并切片
        bool LoadAndSliceModel(const std::string& filepath);

        // 获取给 GPU 准备好的三大件
        const std::vector<MeshletVertex>& GetVertices() const { return m_vertices; }
        const std::vector<uint32_t>& GetMicroIndices() const { return m_microIndices; }
        const std::vector<MeshletDescriptor>& GetMeshlets() const { return m_meshlets; }

    private:
        // 我们最终要塞给 GPU 的三大连续内存块
        std::vector<MeshletVertex> m_vertices;
        std::vector<uint32_t> m_microIndices;
        std::vector<MeshletDescriptor> m_meshlets;
    };

} // namespace Lizeral