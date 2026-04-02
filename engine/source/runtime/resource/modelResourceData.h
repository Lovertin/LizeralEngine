#pragma once

#include "runtime/resource/material/materialAsset.h"
#include "runtime/resource/meshletAssetTypes.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace Lizeral::Resource {

    struct ImportedTextureData {
        std::string sourceTag;
        std::vector<unsigned char> encodedBytes;
    };

    struct ImportedMeshData {
        uint32_t materialIndex {0};
        std::vector<MeshletVertex> vertices;
        std::vector<uint32_t> indices;
    };

    struct ImportedModelData {
        std::string sourcePath;
        std::vector<ImportedMeshData> meshes;
        std::vector<MaterialData> materials;
        std::vector<MaterialAsset> materialAssets;
        std::vector<ImportedTextureData> textures;

        void Clear() {
            sourcePath.clear();
            meshes.clear();
            materials.clear();
            materialAssets.clear();
            textures.clear();
        }
    };

    struct MeshletCookOptions {
        uint32_t globalTextureOffset {0};
        size_t maxVerticesPerMeshlet {64};
        size_t maxTrianglesPerMeshlet {124};
        bool ensureFallbackTextureForTexturelessModel {true};
    };

    struct RuntimeModelData {
        std::vector<MeshletVertex> vertices;
        std::vector<uint32_t> microIndices;
        std::vector<MeshletDescriptor> meshlets;
        std::vector<MeshletBounds> bounds;
        std::vector<std::vector<unsigned char>> allTextures;
        std::vector<MaterialData> materials;
        std::vector<MaterialAsset> materialAssets;

        void Clear() {
            vertices.clear();
            microIndices.clear();
            meshlets.clear();
            bounds.clear();
            allTextures.clear();
            materials.clear();
            materialAssets.clear();
        }

        size_t EstimateCpuBytes() const {
            size_t textureBytes = 0;
            for (const auto& bytes : allTextures) {
                textureBytes += bytes.size();
            }

            return vertices.size() * sizeof(MeshletVertex)
                + microIndices.size() * sizeof(uint32_t)
                + meshlets.size() * sizeof(MeshletDescriptor)
                + bounds.size() * sizeof(MeshletBounds)
                + materials.size() * sizeof(MaterialData)
                + materialAssets.size() * sizeof(MaterialAsset)
                + textureBytes;
        }
    };

} // namespace Lizeral::Resource
