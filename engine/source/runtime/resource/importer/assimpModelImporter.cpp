#include "runtime/resource/importer/assimpModelImporter.h"

#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/scene.h>

#include <cstdlib>
#include <filesystem>
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
            material->pad0 = 0;
            material->pad1 = 0;
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

        Assimp::Importer importer;
        const unsigned int processFlags =
            aiProcess_Triangulate
            | aiProcess_GenNormals
            | aiProcess_FlipUVs
            | aiProcess_PreTransformVertices;

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

        outModel->meshes.reserve(scene->mNumMeshes);
        for (unsigned int meshIndex = 0; meshIndex < scene->mNumMeshes; ++meshIndex) {
            aiMesh* mesh = scene->mMeshes[meshIndex];
            if (!mesh || mesh->mNumVertices == 0 || mesh->mNumFaces == 0) {
                continue;
            }

            ImportedMeshData importedMesh;
            importedMesh.materialIndex = mesh->mMaterialIndex < outModel->materials.size()
                ? mesh->mMaterialIndex
                : 0;

            importedMesh.vertices.reserve(mesh->mNumVertices);
            for (unsigned int vertexIndex = 0; vertexIndex < mesh->mNumVertices; ++vertexIndex) {
                MeshletVertex vertex{};
                vertex.pos[0] = mesh->mVertices[vertexIndex].x;
                vertex.pos[1] = mesh->mVertices[vertexIndex].y;
                vertex.pos[2] = mesh->mVertices[vertexIndex].z;

                if (mesh->HasNormals()) {
                    vertex.normal[0] = mesh->mNormals[vertexIndex].x;
                    vertex.normal[1] = mesh->mNormals[vertexIndex].y;
                    vertex.normal[2] = mesh->mNormals[vertexIndex].z;
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
                outModel->meshes.push_back(std::move(importedMesh));
            }
        }

        if (outModel->meshes.empty()) {
            if (outError) {
                *outError = "Imported scene has no valid triangle mesh.";
            }
            return false;
        }

        std::cout
            << "[ResourceImporter] Imported model: meshes=" << outModel->meshes.size()
            << ", materials=" << outModel->materials.size()
            << ", textures=" << outModel->textures.size()
            << std::endl;

        return true;
    }

} // namespace Lizeral::Resource
