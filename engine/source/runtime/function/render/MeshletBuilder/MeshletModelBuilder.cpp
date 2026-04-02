#include "MeshletModelBuilder.h"

#include "runtime/resource/resourceLoadingUnit.h"

#include <iostream>

namespace Lizeral {

    bool MeshletModelBuilder::LoadAndSliceModel(const std::string& filepath, uint32_t globalTextureOffset) {
        std::cout << "\n[MeshletBuilder] Loading model through resource unit: " << filepath << "..." << std::endl;

        Resource::RuntimeModelData runtimeModelData;
        Resource::ResourceLoadingUnit loadingUnit;

        std::string error;
        if (!loadingUnit.LoadModel(filepath, &runtimeModelData, globalTextureOffset, &error)) {
            std::cerr << "[MeshletBuilder] Resource loading failed: " << error << std::endl;
            return false;
        }

        m_vertices = std::move(runtimeModelData.vertices);
        m_microIndices = std::move(runtimeModelData.microIndices);
        m_meshlets = std::move(runtimeModelData.meshlets);
        m_bounds = std::move(runtimeModelData.bounds);
        m_allTextures = std::move(runtimeModelData.allTextures);
        m_materials = std::move(runtimeModelData.materials);
        m_materialAssets = std::move(runtimeModelData.materialAssets);

        std::cout << "[MeshletBuilder] SUCCESS! Meshlets=" << m_meshlets.size()
                  << ", Vertices=" << m_vertices.size()
                  << ", MicroIndices=" << m_microIndices.size()
                  << ", Materials=" << m_materials.size()
                  << ", Textures=" << m_allTextures.size()
                  << std::endl;
        return true;
    }

} // namespace Lizeral
