#pragma once

#include "runtime/resource/material/materialAsset.h"
#include "runtime/resource/modelResourceData.h"
#include "runtime/resource/meshletAssetTypes.h"

#include <vector>
#include <string>
#include <cstdint>

namespace Lizeral {

    class MeshletModelBuilder {
    public:
        bool LoadAndSliceModel(const std::string& filepath, uint32_t globalTextureOffset = 0);

        const std::vector<MeshletVertex>& GetVertices() const { return m_vertices; }
        const std::vector<uint32_t>& GetMicroIndices() const { return m_microIndices; }
        const std::vector<MeshletDescriptor>& GetMeshlets() const { return m_meshlets; }
        const std::vector<MeshletBounds>& GetBounds() const { return m_bounds; }

        const std::vector<std::vector<unsigned char>>& GetAllTextures() const { return m_allTextures; }
        const std::vector<MaterialData>& GetMaterials() const { return m_materials; }
        const std::vector<Resource::MaterialAsset>& GetMaterialAssets() const { return m_materialAssets; }
        const std::vector<Resource::RuntimeMeshAssetData>& GetMeshAssets() const { return m_meshAssets; }


    private:
        std::vector<MeshletVertex> m_vertices;
        std::vector<uint32_t> m_microIndices;
        std::vector<MeshletDescriptor> m_meshlets;
        std::vector<MeshletBounds> m_bounds;

        std::vector<std::vector<unsigned char>> m_allTextures;
        std::vector<MaterialData> m_materials;
        std::vector<Resource::MaterialAsset> m_materialAssets;
        std::vector<Resource::RuntimeMeshAssetData> m_meshAssets;
    };

} // namespace Lizeral
