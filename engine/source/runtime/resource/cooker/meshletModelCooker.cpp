#include "runtime/resource/cooker/meshletModelCooker.h"

#include <meshoptimizer.h>

#include <iostream>
#include <string>

namespace Lizeral::Resource {

    namespace {

        // Encoded 1x1 RGBA pure white PNG.
        constexpr unsigned char DEFAULT_WHITE_PNG[] = {
            0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
            0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01,
            0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0x15, 0xc4, 0x89, 0x00, 0x00, 0x00,
            0x0d, 0x49, 0x44, 0x41, 0x54, 0x78, 0xda, 0x63, 0xfc, 0xcf, 0xf0, 0x1f,
            0x00, 0x04, 0x01, 0x01, 0x00, 0x61, 0x93, 0x56, 0x22, 0x00, 0x00, 0x00,
            0x00, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
        };

        void PatchTextureIndex(int* textureIndex, uint32_t offset) {
            if (!textureIndex) {
                return;
            }
            if (*textureIndex >= 0) {
                *textureIndex += static_cast<int>(offset);
            }
        }

    } // namespace

    bool MeshletModelCooker::Cook(const ImportedModelData& importedModel, const MeshletCookOptions& options, RuntimeModelData* outModel, std::string* outError) const {
        if (!outModel) {
            if (outError) {
                *outError = "Cook failed: outModel is null.";
            }
            return false;
        }

        outModel->Clear();
        outModel->sourcePath = importedModel.sourcePath;
        outModel->modelName = importedModel.modelName;
        outModel->rootNodeIndex = importedModel.rootNodeIndex;
        outModel->nodes = importedModel.nodes;

        outModel->materials = importedModel.materials;
        outModel->materialAssets = importedModel.materialAssets;
        if (outModel->materialAssets.size() != outModel->materials.size()) {
            outModel->materialAssets.clear();
            outModel->materialAssets.reserve(outModel->materials.size());
            for (size_t i = 0; i < outModel->materials.size(); ++i) {
                outModel->materialAssets.push_back(
                    MakeMaterialAssetFromGpuMaterial(outModel->materials[i], "material_" + std::to_string(i))
                );
            }
        }
        if (outModel->materials.empty()) {
            MaterialData fallbackMaterial{};
            fallbackMaterial.baseColorFactor[0] = 1.0f;
            fallbackMaterial.baseColorFactor[1] = 1.0f;
            fallbackMaterial.baseColorFactor[2] = 1.0f;
            fallbackMaterial.baseColorFactor[3] = 1.0f;
            fallbackMaterial.metallicFactor = 0.0f;
            fallbackMaterial.roughnessFactor = 1.0f;
            fallbackMaterial.albedoTex = -1;
            fallbackMaterial.normalTex = -1;
            fallbackMaterial.ormTex = -1;
            fallbackMaterial.emissiveTex = -1;
            fallbackMaterial.pad0 = 0;
            fallbackMaterial.pad1 = 0;
            outModel->materials.push_back(fallbackMaterial);
            outModel->materialAssets.push_back(MakeMaterialAssetFromGpuMaterial(fallbackMaterial, "material_0"));
        }

        outModel->allTextures.reserve(importedModel.textures.size() + 1);
        for (const auto& texture : importedModel.textures) {
            outModel->allTextures.push_back(texture.encodedBytes);
        }

        if (outModel->allTextures.empty() && options.ensureFallbackTextureForTexturelessModel) {
            outModel->allTextures.emplace_back(DEFAULT_WHITE_PNG, DEFAULT_WHITE_PNG + sizeof(DEFAULT_WHITE_PNG));
        }

        for (auto& material : outModel->materials) {
            PatchTextureIndex(&material.albedoTex, options.globalTextureOffset);
            PatchTextureIndex(&material.normalTex, options.globalTextureOffset);
            PatchTextureIndex(&material.ormTex, options.globalTextureOffset);
            PatchTextureIndex(&material.emissiveTex, options.globalTextureOffset);
        }
        for (auto& materialAsset : outModel->materialAssets) {
            PatchTextureIndex(&materialAsset.textures.albedoTex, options.globalTextureOffset);
            PatchTextureIndex(&materialAsset.textures.normalTex, options.globalTextureOffset);
            PatchTextureIndex(&materialAsset.textures.ormTex, options.globalTextureOffset);
            PatchTextureIndex(&materialAsset.textures.emissiveTex, options.globalTextureOffset);
        }

        for (const auto& mesh : importedModel.meshes) {
            if (mesh.vertices.empty() || mesh.indices.empty()) {
                continue;
            }

            RuntimeMeshAssetData runtimeMeshAsset{};
            runtimeMeshAsset.name = mesh.name;
            runtimeMeshAsset.sourceMeshIndex = mesh.sourceMeshIndex;
            runtimeMeshAsset.nodeIndex = mesh.nodeIndex;
            runtimeMeshAsset.materialIndex = mesh.materialIndex < outModel->materials.size() ? mesh.materialIndex : 0;
            runtimeMeshAsset.meshletOffset = static_cast<uint32_t>(outModel->meshlets.size());
            runtimeMeshAsset.vertexOffset = static_cast<uint32_t>(outModel->vertices.size());
            runtimeMeshAsset.microIndexOffset = static_cast<uint32_t>(outModel->microIndices.size());

            const size_t maxMeshlets = meshopt_buildMeshletsBound(
                mesh.indices.size(),
                options.maxVerticesPerMeshlet,
                options.maxTrianglesPerMeshlet
            );

            std::vector<meshopt_Meshlet> meshlets(maxMeshlets);
            std::vector<unsigned int> meshletVertices(maxMeshlets * options.maxVerticesPerMeshlet);
            std::vector<unsigned char> meshletTriangles(maxMeshlets * options.maxTrianglesPerMeshlet * 3);

            const size_t meshletCount = meshopt_buildMeshlets(
                meshlets.data(),
                meshletVertices.data(),
                meshletTriangles.data(),
                mesh.indices.data(),
                mesh.indices.size(),
                &mesh.vertices[0].pos[0],
                mesh.vertices.size(),
                sizeof(MeshletVertex),
                options.maxVerticesPerMeshlet,
                options.maxTrianglesPerMeshlet,
                0.0f
            );

            meshlets.resize(meshletCount);
            const uint32_t materialIndex = mesh.materialIndex < outModel->materials.size() ? mesh.materialIndex : 0;

            for (const meshopt_Meshlet& rawMeshlet : meshlets) {
                MeshletDescriptor descriptor{};
                descriptor.vertexOffset = static_cast<uint32_t>(outModel->vertices.size());
                descriptor.triangleOffset = static_cast<uint32_t>(outModel->microIndices.size());
                descriptor.vertexCount = rawMeshlet.vertex_count;
                descriptor.triangleCount = rawMeshlet.triangle_count;
                descriptor.materialID = materialIndex;

                const meshopt_Bounds meshletBounds = meshopt_computeMeshletBounds(
                    &meshletVertices[rawMeshlet.vertex_offset],
                    &meshletTriangles[rawMeshlet.triangle_offset],
                    rawMeshlet.triangle_count,
                    &mesh.vertices[0].pos[0],
                    mesh.vertices.size(),
                    sizeof(MeshletVertex)
                );

                MeshletBounds bounds{};
                bounds.center[0] = meshletBounds.center[0];
                bounds.center[1] = meshletBounds.center[1];
                bounds.center[2] = meshletBounds.center[2];
                bounds.radius = meshletBounds.radius;

                outModel->meshlets.push_back(descriptor);
                outModel->bounds.push_back(bounds);

                for (unsigned int vertexIdx = 0; vertexIdx < rawMeshlet.vertex_count; ++vertexIdx) {
                    const unsigned int localVertexIndex = meshletVertices[rawMeshlet.vertex_offset + vertexIdx];
                    outModel->vertices.push_back(mesh.vertices[localVertexIndex]);
                }

                const size_t triangleIndexCount = static_cast<size_t>(rawMeshlet.triangle_count) * 3;
                for (size_t triIdx = 0; triIdx < triangleIndexCount; ++triIdx) {
                    outModel->microIndices.push_back(static_cast<uint32_t>(meshletTriangles[rawMeshlet.triangle_offset + triIdx]));
                }
            }

            runtimeMeshAsset.meshletCount = static_cast<uint32_t>(outModel->meshlets.size()) - runtimeMeshAsset.meshletOffset;
            runtimeMeshAsset.vertexCount = static_cast<uint32_t>(outModel->vertices.size()) - runtimeMeshAsset.vertexOffset;
            runtimeMeshAsset.microIndexCount = static_cast<uint32_t>(outModel->microIndices.size()) - runtimeMeshAsset.microIndexOffset;
            outModel->meshAssets.push_back(std::move(runtimeMeshAsset));
        }

        if (outModel->meshlets.empty()) {
            if (outError) {
                *outError = "Cooked result is empty: no meshlets were generated.";
            }
            return false;
        }

        std::cout
            << "[ResourceCooker] Cooked model: meshlets=" << outModel->meshlets.size()
            << ", vertices=" << outModel->vertices.size()
            << ", microIndices=" << outModel->microIndices.size()
            << ", textures=" << outModel->allTextures.size()
            << ", cpuBytes=" << outModel->EstimateCpuBytes()
            << std::endl;

        return true;
    }

} // namespace Lizeral::Resource
