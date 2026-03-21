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
        uint32_t textureID;    
    };

    struct MaterialData {
        float baseColorFactor[4];
        float metallicFactor;
        float roughnessFactor;
        int albedoTex; 
        int normalTex;
        int ormTex;
        int emissiveTex;
        int pad0, pad1;
    };

    class MeshletModelBuilder {
    public:
        bool LoadAndSliceModel(const std::string& filepath, uint32_t globalTextureOffset = 0);

        const std::vector<MeshletVertex>& GetVertices() const { return m_vertices; }
        const std::vector<uint32_t>& GetMicroIndices() const { return m_microIndices; }
        const std::vector<MeshletDescriptor>& GetMeshlets() const { return m_meshlets; }
        const std::vector<MeshletBounds>& GetBounds() const { return m_bounds; }

        const std::vector<std::vector<unsigned char>>& GetAllTextures() const { return m_allTextures; }
        const std::vector<MaterialData>& GetMaterials() const { return m_materials; }


    private:
        std::vector<MeshletVertex> m_vertices;
        std::vector<uint32_t> m_microIndices;
        std::vector<MeshletDescriptor> m_meshlets;
        std::vector<MeshletBounds> m_bounds;

        std::vector<std::vector<unsigned char>> m_allTextures;
        std::vector<MaterialData> m_materials;
    };

} // namespace Lizeral