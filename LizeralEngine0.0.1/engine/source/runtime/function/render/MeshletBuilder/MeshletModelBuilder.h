#pragma once

#include <vector>
#include <string>
#include <cstdint>

namespace Lizeral {

    struct MeshletVertex {
        float pos[3];
        float normal[3];
        float uv[2];
    };

    struct MeshletBounds {
        float center[3];
        float radius;
    };

    struct MeshletDescriptor {
        uint32_t vertexOffset;    
        uint32_t triangleOffset;  
        uint32_t vertexCount;     
        uint32_t triangleCount;   
        uint32_t textureID;       // ★ 新增：给 Meshlet 发放的贴图身份证！
    };

    struct MaterialData {
        float baseColorFactor[4];
        float metallicFactor;
        float roughnessFactor;
        float padding[2]; // 补齐 16 字节
    };

    class MeshletModelBuilder {
    public:
        bool LoadAndSliceModel(const std::string& filepath);

        const std::vector<MeshletVertex>& GetVertices() const { return m_vertices; }
        const std::vector<uint32_t>& GetMicroIndices() const { return m_microIndices; }
        const std::vector<MeshletDescriptor>& GetMeshlets() const { return m_meshlets; }
        const std::vector<MeshletBounds>& GetBounds() const { return m_bounds; }

        // 新增：获取模型内所有的贴图字节流 (索引对应材质 ID)
        const std::vector<std::vector<unsigned char>>& GetAllTextures() const { return m_allTextures; }
        const std::vector<MaterialData>& GetMaterials() const { return m_materials; }

    private:
        std::vector<MeshletVertex> m_vertices;
        std::vector<uint32_t> m_microIndices;
        std::vector<MeshletDescriptor> m_meshlets;
        std::vector<MeshletBounds> m_bounds;

        // 存储所有材质的贴图内存块
        std::vector<std::vector<unsigned char>> m_allTextures;
        std::vector<MaterialData> m_materials;
    };

} // namespace Lizeral