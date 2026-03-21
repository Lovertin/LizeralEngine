#include "MeshletModelBuilder.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <meshoptimizer.h>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

namespace Lizeral {

    // magic：one 1x1 pure white string
    const unsigned char DEFAULT_WHITE_PNG[] = {
        0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d, 
        0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 
        0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0x15, 0xc4, 0x89, 0x00, 0x00, 0x00, 
        0x0d, 0x49, 0x44, 0x41, 0x54, 0x78, 0xda, 0x63, 0xfc, 0xcf, 0xf0, 0x1f, 
        0x00, 0x04, 0x01, 0x01, 0x00, 0x61, 0x93, 0x56, 0x22, 0x00, 0x00, 0x00, 
        0x00, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
    };

    bool MeshletModelBuilder::LoadAndSliceModel(const std::string& filepath , uint32_t globalTextureOffset) {
        std::cout << "\n[MeshletBuilder] Loading model: " << filepath << "..." << std::endl;

        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile(filepath, 
            aiProcess_Triangulate | aiProcess_GenNormals | 
            aiProcess_FlipUVs | aiProcess_PreTransformVertices
        );

        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "[MeshletBuilder] Assimp ERROR: " << importer.GetErrorString() << std::endl;
            return false;
        }

        m_allTextures.clear();
        m_materials.clear();

        std::unordered_map<int, int> loadedTextures;

        auto LoadTexture = [&](aiMaterial* material, aiTextureType type) -> int {
            aiString texPath;
            if (material->GetTexture(type, 0, &texPath) == AI_SUCCESS) {
                if (texPath.data[0] == '*') {
                    int texIndex = std::stoi(&texPath.data[1]);
                    if (texIndex < scene->mNumTextures) {
                        if (loadedTextures.find(texIndex) != loadedTextures.end()) {
                            return loadedTextures[texIndex];
                        }
                        
                        aiTexture* aiTex = scene->mTextures[texIndex];
                        if (aiTex->mHeight == 0) { //  PNG / JPG
                            int byteLength = aiTex->mWidth;
                            const unsigned char* rawData = reinterpret_cast<const unsigned char*>(aiTex->pcData);
                            
                            m_allTextures.push_back(std::vector<unsigned char>(rawData, rawData + byteLength));
                            int newLocalIndex = static_cast<int>(m_allTextures.size() - 1);
                            
                            loadedTextures[texIndex] = newLocalIndex; //register
                            return newLocalIndex;
                        }
                    }
                }
            }
            return -1; // find null
        };

        for (unsigned int i = 0; i < scene->mNumMaterials; i++) {
            aiMaterial* material = scene->mMaterials[i];
            MaterialData matData{};

            // set pbr default
            matData.baseColorFactor[0] = 1.0f; matData.baseColorFactor[1] = 1.0f; 
            matData.baseColorFactor[2] = 1.0f; matData.baseColorFactor[3] = 1.0f;
            matData.metallicFactor = 1.0f;  // glTF standard
            matData.roughnessFactor = 1.0f;

            aiColor4D baseColor;
            if (AI_SUCCESS == material->Get(AI_MATKEY_BASE_COLOR, baseColor)) {
                matData.baseColorFactor[0] = baseColor.r; matData.baseColorFactor[1] = baseColor.g;
                matData.baseColorFactor[2] = baseColor.b; matData.baseColorFactor[3] = baseColor.a;
            }
            material->Get(AI_MATKEY_METALLIC_FACTOR, matData.metallicFactor);
            material->Get(AI_MATKEY_ROUGHNESS_FACTOR, matData.roughnessFactor);

            // Albedo 
            int albedoLocal = LoadTexture(material, aiTextureType_BASE_COLOR);
            if (albedoLocal == -1) albedoLocal = LoadTexture(material, aiTextureType_DIFFUSE);
            matData.albedoTex = (albedoLocal != -1) ? (globalTextureOffset + albedoLocal) : -1;

            // Normal
            int normalLocal = LoadTexture(material, aiTextureType_NORMALS);
            matData.normalTex = (normalLocal != -1) ? (globalTextureOffset + normalLocal) : -1;

            //pick ORM (AO/Roughness/Metallic) 
            int ormLocal = LoadTexture(material, aiTextureType_UNKNOWN); 
            if (ormLocal == -1) ormLocal = LoadTexture(material, aiTextureType_METALNESS);
            if (ormLocal == -1) ormLocal = LoadTexture(material, aiTextureType_DIFFUSE_ROUGHNESS);
            matData.ormTex = (ormLocal != -1) ? (globalTextureOffset + ormLocal) : -1;

            // pick Emissive 
            int emissiveLocal = LoadTexture(material, aiTextureType_EMISSION_COLOR);
            if (emissiveLocal == -1) emissiveLocal = LoadTexture(material, aiTextureType_EMISSIVE);
            matData.emissiveTex = (emissiveLocal != -1) ? (globalTextureOffset + emissiveLocal) : -1;

            matData.pad0 = 0; matData.pad1 = 0;
            m_materials.push_back(matData);
        }
        
        std::cout << "[MeshletBuilder] Extracted " << m_allTextures.size() << " UNIQUE textures for " << m_materials.size() << " materials." << std::endl;

        m_vertices.clear();
        m_microIndices.clear();
        m_meshlets.clear();
        m_bounds.clear();

        for (unsigned int i = 0; i < scene->mNumMeshes; i++) {
            aiMesh* mesh = scene->mMeshes[i];
            uint32_t materialID = mesh->mMaterialIndex;

            std::vector<MeshletVertex> localVertices;
            std::vector<uint32_t> localIndices;

            for (unsigned int j = 0; j < mesh->mNumVertices; j++) {
                MeshletVertex v{};
                v.pos[0] = mesh->mVertices[j].x;
                v.pos[1] = mesh->mVertices[j].y;
                v.pos[2] = mesh->mVertices[j].z;
                if (mesh->HasNormals()) {
                    v.normal[0] = mesh->mNormals[j].x;
                    v.normal[1] = mesh->mNormals[j].y;
                    v.normal[2] = mesh->mNormals[j].z;
                }
                if (mesh->mTextureCoords[0]) {
                    v.uv[0] = mesh->mTextureCoords[0][j].x;
                    v.uv[1] = mesh->mTextureCoords[0][j].y;
                }
                localVertices.push_back(v);
            }

            for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
                aiFace face = mesh->mFaces[j];
                for (unsigned int k = 0; k < face.mNumIndices; k++) {
                    localIndices.push_back(face.mIndices[k]);
                }
            }

            // slicing this submeshes
            const size_t max_vertices = 64;
            const size_t max_triangles = 124; 
            size_t max_meshlets = meshopt_buildMeshletsBound(localIndices.size(), max_vertices, max_triangles);
            
            std::vector<meshopt_Meshlet> localMeshlets(max_meshlets);
            std::vector<unsigned int> meshletVertices(max_meshlets * max_vertices);
            std::vector<unsigned char> meshletTriangles(max_meshlets * max_triangles * 3);

            size_t meshlet_count = meshopt_buildMeshlets(
                localMeshlets.data(), meshletVertices.data(), meshletTriangles.data(),
                localIndices.data(), localIndices.size(), 
                &localVertices[0].pos[0], localVertices.size(), sizeof(MeshletVertex), 
                max_vertices, max_triangles, 0.0f);

            localMeshlets.resize(meshlet_count);

            uint32_t currentGlobalVertexCount = static_cast<uint32_t>(m_vertices.size());
            uint32_t currentGlobalTriangleCount = static_cast<uint32_t>(m_microIndices.size());

            for (const auto& m : localMeshlets) {
                MeshletDescriptor descriptor{};
                descriptor.vertexOffset = currentGlobalVertexCount + static_cast<uint32_t>(m_vertices.size() - currentGlobalVertexCount);
                descriptor.triangleOffset = currentGlobalTriangleCount + static_cast<uint32_t>(m_microIndices.size() - currentGlobalTriangleCount);
                descriptor.vertexCount = m.vertex_count;
                descriptor.triangleCount = m.triangle_count;
                descriptor.textureID = globalTextureOffset + materialID;
                m_meshlets.push_back(descriptor);

                meshopt_Bounds bounds = meshopt_computeMeshletBounds(
                    &meshletVertices[m.vertex_offset],     
                    &meshletTriangles[m.triangle_offset],  
                    m.triangle_count,                      
                    &localVertices[0].pos[0],             
                    localVertices.size(),                 
                    sizeof(MeshletVertex)                  
                );

                MeshletBounds myBounds;
                myBounds.center[0] = bounds.center[0];
                myBounds.center[1] = bounds.center[1];
                myBounds.center[2] = bounds.center[2];
                myBounds.radius    = bounds.radius;
                m_bounds.push_back(myBounds);

                for (unsigned int vIdx = 0; vIdx < m.vertex_count; ++vIdx) {
                    m_vertices.push_back(localVertices[meshletVertices[m.vertex_offset + vIdx]]);
                }

                for (unsigned int tIdx = 0; tIdx < m.triangle_count * 3; ++tIdx) {
                    m_microIndices.push_back(static_cast<uint32_t>(meshletTriangles[m.triangle_offset + tIdx]));
                }
            }
        }

        std::cout << "[MeshletBuilder] SUCCESS! Model sliced into " << m_meshlets.size() << " highly optimized Meshlets!" << std::endl;
        std::cout << "[MeshletBuilder] Final buffers: Vertices(" << m_vertices.size() << "), MicroIndices(" << m_microIndices.size() << "), Bounds(" << m_bounds.size() << ")" << std::endl;

        return true;
    }

} // namespace Lizeral