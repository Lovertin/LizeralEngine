#pragma once

#include "runtime/resource/material/materialAsset.h"
#include "runtime/resource/meshletAssetTypes.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace Lizeral::Resource {

    struct NodeTransformData {
        float translation[3] {0.0f, 0.0f, 0.0f};
        float rotation[4] {0.0f, 0.0f, 0.0f, 1.0f};
        float scale[3] {1.0f, 1.0f, 1.0f};
    };

    struct ImportedTextureData {
        std::string name;
        std::string sourceTag;
        std::vector<unsigned char> encodedBytes;
    };

    struct ImportedMeshData {
        std::string name;
        uint32_t sourceMeshIndex {0};
        uint32_t nodeIndex {0};
        uint32_t materialIndex {0};
        std::vector<MeshletVertex> vertices;
        std::vector<uint32_t> indices;
    };

    struct ModelNodeData {
        std::string name;
        int32_t parentIndex {-1};
        NodeTransformData localTransform {};
        std::vector<uint32_t> childNodeIndices;
        std::vector<uint32_t> meshIndices;
    };

    struct RuntimeMeshAssetData {
        std::string name;
        uint32_t sourceMeshIndex {0};
        uint32_t nodeIndex {0};
        uint32_t materialIndex {0};
        uint32_t meshletOffset {0};
        uint32_t meshletCount {0};
        uint32_t vertexOffset {0};
        uint32_t vertexCount {0};
        uint32_t microIndexOffset {0};
        uint32_t microIndexCount {0};
    };

    struct ImportedModelData {
        std::string sourcePath;
        std::string modelName;
        int32_t rootNodeIndex {-1};
        std::vector<ModelNodeData> nodes;
        std::vector<ImportedMeshData> meshes;
        std::vector<MaterialData> materials;
        std::vector<MaterialAsset> materialAssets;
        std::vector<ImportedTextureData> textures;

        void Clear() {
            sourcePath.clear();
            modelName.clear();
            rootNodeIndex = -1;
            nodes.clear();
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
        std::string sourcePath;
        std::string modelName;
        int32_t rootNodeIndex {-1};
        std::vector<ModelNodeData> nodes;
        std::vector<RuntimeMeshAssetData> meshAssets;
        std::vector<MeshletVertex> vertices;
        std::vector<uint32_t> microIndices;
        std::vector<MeshletDescriptor> meshlets;
        std::vector<MeshletBounds> bounds;
        std::vector<std::vector<unsigned char>> allTextures;
        std::vector<MaterialData> materials;
        std::vector<MaterialAsset> materialAssets;

        void Clear() {
            sourcePath.clear();
            modelName.clear();
            rootNodeIndex = -1;
            nodes.clear();
            meshAssets.clear();
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
                + nodes.size() * sizeof(ModelNodeData)
                + meshAssets.size() * sizeof(RuntimeMeshAssetData)
                + materials.size() * sizeof(MaterialData)
                + materialAssets.size() * sizeof(MaterialAsset)
                + textureBytes;
        }
    };

} // namespace Lizeral::Resource
