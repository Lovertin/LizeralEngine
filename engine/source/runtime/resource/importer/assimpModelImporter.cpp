#include "runtime/resource/importer/assimpModelImporter.h"

#include <assimp/Importer.hpp>
#include <assimp/GltfMaterial.h>
#include <assimp/postprocess.h>
#include <assimp/scene.h>

#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <unordered_map>

namespace Lizeral::Resource {

    namespace {

        bool ReadBinaryFile(const std::filesystem::path& path, std::vector<unsigned char>* outBytes) {
            if (!outBytes) {
                return false;
            }

            std::ifstream file(path, std::ios::binary | std::ios::ate);
            if (!file.is_open()) {
                return false;
            }

            const std::streamsize size = file.tellg();
            if (size < 0) {
                return false;
            }

            outBytes->resize(static_cast<size_t>(size));
            file.seekg(0, std::ios::beg);
            if (size > 0) {
                file.read(reinterpret_cast<char*>(outBytes->data()), size);
            }

            return file.good() || file.eof();
        }

        std::string BuildTextureKey(const std::filesystem::path& path) {
            std::error_code ec;
            const auto canonical = std::filesystem::weakly_canonical(path, ec);
            if (!ec) {
                return canonical.generic_string();
            }
            return path.lexically_normal().generic_string();
        }

        void SetDefaultMaterial(MaterialData* material) {
            if (!material) {
                return;
            }

            material->baseColorFactor[0] = 1.0f;
            material->baseColorFactor[1] = 1.0f;
            material->baseColorFactor[2] = 1.0f;
            material->baseColorFactor[3] = 1.0f;
            // Most DCC/PBR assets are non-metal by default.
            material->metallicFactor = 0.0f;
            material->roughnessFactor = 1.0f;

            material->albedoTex = -1;
            material->normalTex = -1;
            material->ormTex = -1;
            material->emissiveTex = -1;
            material->alphaMode = static_cast<int>(MaterialAlphaMode::Opaque);
            material->alphaCutoff = 0.5f;
            material->pad0 = 0;
            material->pad1 = 0;
        }

        void FillNodeTransform(const aiMatrix4x4& transform, NodeTransformData* outTransform) {
            if (!outTransform) {
                return;
            }

            aiVector3D scaling;
            aiQuaternion rotation;
            aiVector3D translation;
            transform.Decompose(scaling, rotation, translation);

            outTransform->translation[0] = translation.x;
            outTransform->translation[1] = translation.y;
            outTransform->translation[2] = translation.z;

            outTransform->rotation[0] = rotation.x;
            outTransform->rotation[1] = rotation.y;
            outTransform->rotation[2] = rotation.z;
            outTransform->rotation[3] = rotation.w;

            outTransform->scale[0] = scaling.x;
            outTransform->scale[1] = scaling.y;
            outTransform->scale[2] = scaling.z;
        }

        std::string BuildMeshDisplayName(const aiMesh* mesh, const aiNode* ownerNode, unsigned int sourceMeshIndex) {
            if (mesh && mesh->mName.length > 0) {
                return mesh->mName.C_Str();
            }

            if (ownerNode && ownerNode->mName.length > 0) {
                return std::string(ownerNode->mName.C_Str()) + "_mesh_" + std::to_string(sourceMeshIndex);
            }

            return "mesh_" + std::to_string(sourceMeshIndex);
        }

    } // namespace

    bool AssimpModelImporter::ImportModel(const std::string& sourcePath, ImportedModelData* outModel, std::string* outError) const {
        if (!outModel) {
            if (outError) {
                *outError = "ImportModel failed: outModel is null.";
            }
            return false;
        }

        outModel->Clear();
        outModel->sourcePath = sourcePath;
        {
            std::filesystem::path sourceFsPath(sourcePath);
            outModel->modelName = sourceFsPath.stem().string();
            if (outModel->modelName.empty()) {
                outModel->modelName = "model";
            }
        }

        Assimp::Importer importer;
        const unsigned int processFlags =
            aiProcess_Triangulate
            | aiProcess_GenNormals
            | aiProcess_FlipUVs;

        const aiScene* scene = importer.ReadFile(sourcePath, processFlags);
        if (!scene || (scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE) || !scene->mRootNode) {
            if (outError) {
                *outError = importer.GetErrorString();
            }
            return false;
        }

        std::filesystem::path baseDir;
        {
            std::error_code ec;
            auto absPath = std::filesystem::absolute(std::filesystem::path(sourcePath), ec);
            if (ec) {
                absPath = std::filesystem::path(sourcePath);
            }
            baseDir = absPath.parent_path();
        }

        std::unordered_map<int, int> embeddedTextureToLocalIndex;
        std::unordered_map<std::string, int> externalTextureToLocalIndex;

        auto loadTextureRef = [&](aiMaterial* material, aiTextureType type) -> int {
            if (!material || material->GetTextureCount(type) == 0) {
                return -1;
            }

            aiString texturePath;
            if (material->GetTexture(type, 0, &texturePath) != AI_SUCCESS) {
                return -1;
            }

            if (texturePath.length == 0) {
                return -1;
            }

            const std::string texturePathStr(texturePath.C_Str());

            if (!texturePathStr.empty() && texturePathStr[0] == '*') {
                const int embeddedIndex = std::atoi(texturePathStr.c_str() + 1);
                if (embeddedIndex < 0 || static_cast<unsigned int>(embeddedIndex) >= scene->mNumTextures) {
                    return -1;
                }

                auto it = embeddedTextureToLocalIndex.find(embeddedIndex);
                if (it != embeddedTextureToLocalIndex.end()) {
                    return it->second;
                }

                const aiTexture* embeddedTexture = scene->mTextures[embeddedIndex];
                if (!embeddedTexture || embeddedTexture->mHeight != 0) {
                    // mHeight != 0 means raw pixel block instead of encoded image stream.
                    return -1;
                }

                const auto* rawData = reinterpret_cast<const unsigned char*>(embeddedTexture->pcData);
                const int byteCount = embeddedTexture->mWidth;
                if (!rawData || byteCount <= 0) {
                    return -1;
                }

                ImportedTextureData importedTexture;
                importedTexture.name = "embedded_" + std::to_string(embeddedIndex);
                importedTexture.sourceTag = "embedded:" + std::to_string(embeddedIndex);
                importedTexture.encodedBytes.assign(rawData, rawData + byteCount);

                outModel->textures.push_back(std::move(importedTexture));
                const int localTextureIndex = static_cast<int>(outModel->textures.size() - 1);
                embeddedTextureToLocalIndex[embeddedIndex] = localTextureIndex;
                return localTextureIndex;
            }

            std::filesystem::path candidate = std::filesystem::path(texturePathStr);
            if (candidate.is_relative()) {
                candidate = baseDir / candidate;
            }

            if (!std::filesystem::exists(candidate)) {
                std::filesystem::path fallbackCandidate = std::filesystem::path(texturePathStr);
                if (std::filesystem::exists(fallbackCandidate)) {
                    candidate = fallbackCandidate;
                } else {
                    return -1;
                }
            }

            const std::string key = BuildTextureKey(candidate);
            auto existingTexture = externalTextureToLocalIndex.find(key);
            if (existingTexture != externalTextureToLocalIndex.end()) {
                return existingTexture->second;
            }

            ImportedTextureData importedTexture;
            importedTexture.name = candidate.stem().string();
            if (importedTexture.name.empty()) {
                importedTexture.name = "texture_" + std::to_string(outModel->textures.size());
            }
            importedTexture.sourceTag = key;
            if (!ReadBinaryFile(candidate, &importedTexture.encodedBytes)) {
                return -1;
            }

            outModel->textures.push_back(std::move(importedTexture));
            const int localTextureIndex = static_cast<int>(outModel->textures.size() - 1);
            externalTextureToLocalIndex[key] = localTextureIndex;
            return localTextureIndex;
        };

        outModel->materials.reserve(scene->mNumMaterials > 0 ? scene->mNumMaterials : 1);
        outModel->materialAssets.reserve(scene->mNumMaterials > 0 ? scene->mNumMaterials : 1);
        for (unsigned int materialIndex = 0; materialIndex < scene->mNumMaterials; ++materialIndex) {
            aiMaterial* material = scene->mMaterials[materialIndex];
            MaterialData importedMaterial{};
            SetDefaultMaterial(&importedMaterial);

            aiColor4D baseColor;
            if (AI_SUCCESS == material->Get(AI_MATKEY_BASE_COLOR, baseColor)) {
                importedMaterial.baseColorFactor[0] = baseColor.r;
                importedMaterial.baseColorFactor[1] = baseColor.g;
                importedMaterial.baseColorFactor[2] = baseColor.b;
                importedMaterial.baseColorFactor[3] = baseColor.a;
            } else if (AI_SUCCESS == material->Get(AI_MATKEY_COLOR_DIFFUSE, baseColor)) {
                importedMaterial.baseColorFactor[0] = baseColor.r;
                importedMaterial.baseColorFactor[1] = baseColor.g;
                importedMaterial.baseColorFactor[2] = baseColor.b;
                importedMaterial.baseColorFactor[3] = 1.0f;
            }

            material->Get(AI_MATKEY_METALLIC_FACTOR, importedMaterial.metallicFactor);
            material->Get(AI_MATKEY_ROUGHNESS_FACTOR, importedMaterial.roughnessFactor);

            aiString alphaModeStr;
            if (AI_SUCCESS == material->Get(AI_MATKEY_GLTF_ALPHAMODE, alphaModeStr) && alphaModeStr.length > 0) {
                const std::string mode = alphaModeStr.C_Str();
                if (mode == "MASK") {
                    importedMaterial.alphaMode = static_cast<int>(MaterialAlphaMode::Mask);
                } else if (mode == "BLEND") {
                    importedMaterial.alphaMode = static_cast<int>(MaterialAlphaMode::Blend);
                }
            }

            material->Get(AI_MATKEY_GLTF_ALPHACUTOFF, importedMaterial.alphaCutoff);

            float opacity = 1.0f;
            if (AI_SUCCESS == material->Get(AI_MATKEY_OPACITY, opacity)) {
                importedMaterial.baseColorFactor[3] *= opacity;
                if (opacity < 0.999f && importedMaterial.alphaMode == static_cast<int>(MaterialAlphaMode::Opaque)) {
                    importedMaterial.alphaMode = static_cast<int>(MaterialAlphaMode::Blend);
                }
            }

            importedMaterial.albedoTex = loadTextureRef(material, aiTextureType_BASE_COLOR);
            if (importedMaterial.albedoTex < 0) {
                importedMaterial.albedoTex = loadTextureRef(material, aiTextureType_DIFFUSE);
            }

            importedMaterial.normalTex = loadTextureRef(material, aiTextureType_NORMALS);
            if (importedMaterial.normalTex < 0) {
                importedMaterial.normalTex = loadTextureRef(material, aiTextureType_NORMAL_CAMERA);
            }

            importedMaterial.ormTex = loadTextureRef(material, aiTextureType_METALNESS);
            if (importedMaterial.ormTex < 0) {
                importedMaterial.ormTex = loadTextureRef(material, aiTextureType_DIFFUSE_ROUGHNESS);
            }
            if (importedMaterial.ormTex < 0) {
                importedMaterial.ormTex = loadTextureRef(material, aiTextureType_AMBIENT_OCCLUSION);
            }
            if (importedMaterial.ormTex < 0) {
                importedMaterial.ormTex = loadTextureRef(material, aiTextureType_UNKNOWN);
            }

            importedMaterial.emissiveTex = loadTextureRef(material, aiTextureType_EMISSION_COLOR);
            if (importedMaterial.emissiveTex < 0) {
                importedMaterial.emissiveTex = loadTextureRef(material, aiTextureType_EMISSIVE);
            }

            outModel->materials.push_back(importedMaterial);

            aiString materialName;
            std::string materialNameStr;
            if (AI_SUCCESS == material->Get(AI_MATKEY_NAME, materialName) && materialName.length > 0) {
                materialNameStr = materialName.C_Str();
            } else {
                materialNameStr = "material_" + std::to_string(materialIndex);
            }
            outModel->materialAssets.push_back(MakeMaterialAssetFromGpuMaterial(importedMaterial, materialNameStr));
        }

        if (outModel->materials.empty()) {
            MaterialData fallbackMaterial{};
            SetDefaultMaterial(&fallbackMaterial);
            outModel->materials.push_back(fallbackMaterial);
            outModel->materialAssets.push_back(MakeMaterialAssetFromGpuMaterial(fallbackMaterial, "material_0"));
        }

        outModel->nodes.reserve(scene->mNumMeshes + 1);
        outModel->meshes.reserve(scene->mNumMeshes);

        std::function<uint32_t(aiNode*, int32_t, const aiMatrix4x4&)> processNode =
            [&](aiNode* node, int32_t parentNodeIndex, const aiMatrix4x4& parentTransform) -> uint32_t {
                const uint32_t nodeIndex = static_cast<uint32_t>(outModel->nodes.size());
                outModel->nodes.emplace_back();

                ModelNodeData& nodeData = outModel->nodes.back();
                nodeData.name = (node && node->mName.length > 0) ? node->mName.C_Str() : ("node_" + std::to_string(nodeIndex));
                nodeData.parentIndex = parentNodeIndex;
                FillNodeTransform(node->mTransformation, &nodeData.localTransform);

                if (parentNodeIndex >= 0) {
                    outModel->nodes[static_cast<size_t>(parentNodeIndex)].childNodeIndices.push_back(nodeIndex);
                }

                const aiMatrix4x4 accumulatedTransform = parentTransform * node->mTransformation;
                aiMatrix3x3 normalTransform(accumulatedTransform);
                normalTransform.Inverse();
                normalTransform.Transpose();

                for (unsigned int meshRefIndex = 0; meshRefIndex < node->mNumMeshes; ++meshRefIndex) {
                    const unsigned int sourceMeshIndex = node->mMeshes[meshRefIndex];
                    aiMesh* mesh = scene->mMeshes[sourceMeshIndex];
                    if (!mesh || mesh->mNumVertices == 0 || mesh->mNumFaces == 0) {
                        continue;
                    }

                    ImportedMeshData importedMesh;
                    importedMesh.name = BuildMeshDisplayName(mesh, node, sourceMeshIndex);
                    importedMesh.sourceMeshIndex = sourceMeshIndex;
                    importedMesh.nodeIndex = nodeIndex;
                    importedMesh.materialIndex = mesh->mMaterialIndex < outModel->materials.size()
                        ? mesh->mMaterialIndex
                        : 0;

                    importedMesh.vertices.reserve(mesh->mNumVertices);
                    for (unsigned int vertexIndex = 0; vertexIndex < mesh->mNumVertices; ++vertexIndex) {
                        MeshletVertex vertex{};

                        const aiVector3D transformedPosition = accumulatedTransform * mesh->mVertices[vertexIndex];
                        vertex.pos[0] = transformedPosition.x;
                        vertex.pos[1] = transformedPosition.y;
                        vertex.pos[2] = transformedPosition.z;

                        if (mesh->HasNormals()) {
                            aiVector3D transformedNormal = normalTransform * mesh->mNormals[vertexIndex];
                            transformedNormal.Normalize();
                            vertex.normal[0] = transformedNormal.x;
                            vertex.normal[1] = transformedNormal.y;
                            vertex.normal[2] = transformedNormal.z;
                        }

                        if (mesh->mTextureCoords[0]) {
                            vertex.uv[0] = mesh->mTextureCoords[0][vertexIndex].x;
                            vertex.uv[1] = mesh->mTextureCoords[0][vertexIndex].y;
                        }

                        importedMesh.vertices.push_back(vertex);
                    }

                    importedMesh.indices.reserve(static_cast<size_t>(mesh->mNumFaces) * 3);
                    for (unsigned int faceIndex = 0; faceIndex < mesh->mNumFaces; ++faceIndex) {
                        const aiFace& face = mesh->mFaces[faceIndex];
                        if (face.mNumIndices != 3) {
                            continue;
                        }

                        importedMesh.indices.push_back(face.mIndices[0]);
                        importedMesh.indices.push_back(face.mIndices[1]);
                        importedMesh.indices.push_back(face.mIndices[2]);
                    }

                    if (!importedMesh.indices.empty()) {
                        const uint32_t meshIndex = static_cast<uint32_t>(outModel->meshes.size());
                        outModel->meshes.push_back(std::move(importedMesh));
                        nodeData.meshIndices.push_back(meshIndex);
                    }
                }

                for (unsigned int childIndex = 0; childIndex < node->mNumChildren; ++childIndex) {
                    processNode(node->mChildren[childIndex], static_cast<int32_t>(nodeIndex), accumulatedTransform);
                }

                return nodeIndex;
            };

        aiMatrix4x4 identityTransform;
        outModel->rootNodeIndex = static_cast<int32_t>(processNode(scene->mRootNode, -1, identityTransform));

        if (outModel->meshes.empty()) {
            if (outError) {
                *outError = "Imported scene has no valid triangle mesh.";
            }
            return false;
        }

        std::cout
            << "[ResourceImporter] Imported model '" << outModel->modelName << "'"
            << ": nodes=" << outModel->nodes.size()
            << ", meshes=" << outModel->meshes.size()
            << ", materials=" << outModel->materials.size()
            << ", textures=" << outModel->textures.size()
            << std::endl;

        return true;
    }

} // namespace Lizeral::Resource
