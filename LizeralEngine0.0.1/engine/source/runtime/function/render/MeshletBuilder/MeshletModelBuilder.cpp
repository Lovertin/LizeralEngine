#include "MeshletModelBuilder.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <meshoptimizer.h>
#include <iostream>
#include <stdexcept>

namespace Lizeral {

    // 魔法：一段 1x1 纯白 PNG 图片的十六进制字节流！
    // 用途：遇到没有贴图的材质时，把它当做默认贴图塞给 GPU，防止管线崩溃。
    const unsigned char DEFAULT_WHITE_PNG[] = {
        0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d, 
        0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 
        0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0x15, 0xc4, 0x89, 0x00, 0x00, 0x00, 
        0x0d, 0x49, 0x44, 0x41, 0x54, 0x78, 0xda, 0x63, 0xfc, 0xcf, 0xf0, 0x1f, 
        0x00, 0x04, 0x01, 0x01, 0x00, 0x61, 0x93, 0x56, 0x22, 0x00, 0x00, 0x00, 
        0x00, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
    };

    bool MeshletModelBuilder::LoadAndSliceModel(const std::string& filepath) {
        std::cout << "\n[MeshletBuilder] Loading model: " << filepath << "..." << std::endl;

        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile(filepath, 
            aiProcess_Triangulate | aiProcess_GenSmoothNormals | 
            aiProcess_FlipUVs | aiProcess_PreTransformVertices
        );

        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "[MeshletBuilder] Assimp ERROR: " << importer.GetErrorString() << std::endl;
            return false;
        }

        // ====================================================================
        // ★ 阶段 0：提取所有的材质贴图，构建材质 ID 映射库
        // ====================================================================
        m_allTextures.clear();
        m_allTextures.resize(scene->mNumMaterials);

        for (unsigned int i = 0; i < scene->mNumMaterials; i++) {
            aiMaterial* material = scene->mMaterials[i];
            aiString texPath;
            bool foundTexture = false;

            if (material->GetTexture(aiTextureType_BASE_COLOR, 0, &texPath) == AI_SUCCESS ||
                material->GetTexture(aiTextureType_DIFFUSE, 0, &texPath) == AI_SUCCESS) {
                
                if (texPath.data[0] == '*') {
                    int texIndex = std::stoi(&texPath.data[1]);
                    if (texIndex < scene->mNumTextures) {
                        aiTexture* aiTex = scene->mTextures[texIndex];
                        if (aiTex->mHeight == 0) { // 代表是压缩的二进制流
                            int byteLength = aiTex->mWidth;
                            const unsigned char* rawData = reinterpret_cast<const unsigned char*>(aiTex->pcData);
                            m_allTextures[i].assign(rawData, rawData + byteLength);
                            foundTexture = true;
                        }
                    }
                }
            }

            if (!foundTexture) {
                // 没有贴图的材质？安排兜底纯白贴图！
                m_allTextures[i].assign(std::begin(DEFAULT_WHITE_PNG), std::end(DEFAULT_WHITE_PNG));
            }

            MaterialData matData{};
            // 默认值：纯白、粗糙、非金属
            matData.baseColorFactor[0] = 1.0f; matData.baseColorFactor[1] = 1.0f; 
            matData.baseColorFactor[2] = 1.0f; matData.baseColorFactor[3] = 1.0f;
            matData.metallicFactor = 0.0f;
            matData.roughnessFactor = 1.0f;

            aiColor4D baseColor;
            if (AI_SUCCESS == material->Get(AI_MATKEY_BASE_COLOR, baseColor)) {
                matData.baseColorFactor[0] = baseColor.r;
                matData.baseColorFactor[1] = baseColor.g;
                matData.baseColorFactor[2] = baseColor.b;
                matData.baseColorFactor[3] = baseColor.a;
            }
            // 尝试读取金属度和粗糙度（如果没有，就保持默认）
            material->Get(AI_MATKEY_METALLIC_FACTOR, matData.metallicFactor);
            material->Get(AI_MATKEY_ROUGHNESS_FACTOR, matData.roughnessFactor);

            m_materials.push_back(matData);
        }
        std::cout << "[MeshletBuilder] Extracted " << m_allTextures.size() << " materials/textures." << std::endl;

        // ====================================================================
        // ★ 阶段 1 & 2 & 3：按子网格 (Mesh) 独立切分，死守材质边界！
        // ====================================================================
        m_vertices.clear();
        m_microIndices.clear();
        m_meshlets.clear();
        m_bounds.clear();

        for (unsigned int i = 0; i < scene->mNumMeshes; i++) {
            aiMesh* mesh = scene->mMeshes[i];
            uint32_t materialID = mesh->mMaterialIndex; // ★ 这个子网格的身份证

            std::vector<MeshletVertex> localVertices;
            std::vector<uint32_t> localIndices;

            // 1. 提取当前局部网格的顶点
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

            // 2. 提取当前局部网格的索引
            for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
                aiFace face = mesh->mFaces[j];
                for (unsigned int k = 0; k < face.mNumIndices; k++) {
                    localIndices.push_back(face.mIndices[k]);
                }
            }

            // 3. 呼叫神库！单独切分这个子网格！
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

            // 4. 将局部切分结果，拼接到全局的巨型 Buffer 里！
            // 记住目前的全局兵力（Offset）
            uint32_t currentGlobalVertexCount = static_cast<uint32_t>(m_vertices.size());
            uint32_t currentGlobalTriangleCount = static_cast<uint32_t>(m_microIndices.size());

            for (const auto& m : localMeshlets) {
                MeshletDescriptor descriptor{};
                // ★ 极其关键的寻址魔法：全局已有的总数 + 局部数组在生成时的累加数
                descriptor.vertexOffset = currentGlobalVertexCount + static_cast<uint32_t>(m_vertices.size() - currentGlobalVertexCount);
                descriptor.triangleOffset = currentGlobalTriangleCount + static_cast<uint32_t>(m_microIndices.size() - currentGlobalTriangleCount);
                descriptor.vertexCount = m.vertex_count;
                descriptor.triangleCount = m.triangle_count;
                descriptor.textureID = materialID; // ★ 给这个 Meshlet 发放终身绑定的贴图 ID！
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

                // 将局部的物理顶点复制到全局
                for (unsigned int vIdx = 0; vIdx < m.vertex_count; ++vIdx) {
                    m_vertices.push_back(localVertices[meshletVertices[m.vertex_offset + vIdx]]);
                }

                // 将局部的微型索引复制到全局
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