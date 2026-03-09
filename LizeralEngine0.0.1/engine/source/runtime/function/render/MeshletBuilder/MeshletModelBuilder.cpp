#include "MeshletModelBuilder.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <meshoptimizer.h>
#include <iostream>
#include <stdexcept>

namespace Lizeral {

    bool MeshletModelBuilder::LoadAndSliceModel(const std::string& filepath) {
        std::cout << "\n[MeshletBuilder] Loading model: " << filepath << "..." << std::endl;

        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile(filepath, 
            aiProcess_Triangulate | 
            aiProcess_GenSmoothNormals | 
            aiProcess_FlipUVs | 
            aiProcess_PreTransformVertices // ★ 关键魔法：把所有子网格的 Transform 提前算好并拍扁！
        );

        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "[MeshletBuilder] Assimp ERROR: " << importer.GetErrorString() << std::endl;
            return false;
        }

        // ====================================================================
        // 阶段 1：将整个场景的所有子 Mesh“拍扁”合并成一个巨型数组
        // ====================================================================
        std::vector<MeshletVertex> globalVertices;
        std::vector<uint32_t> globalIndices;

        for (unsigned int i = 0; i < scene->mNumMeshes; i++) {
            aiMesh* mesh = scene->mMeshes[i];
            uint32_t vertexOffset = static_cast<uint32_t>(globalVertices.size());

            // 提取顶点
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
                globalVertices.push_back(v);
            }

            // 提取索引 (加上偏移量)
            for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
                aiFace face = mesh->mFaces[j];
                for (unsigned int k = 0; k < face.mNumIndices; k++) {
                    globalIndices.push_back(face.mIndices[k] + vertexOffset);
                }
            }
        }

        std::cout << "[MeshletBuilder] Flattened model into " << globalVertices.size() << " vertices and " 
                  << (globalIndices.size() / 3) << " triangles." << std::endl;

        // ====================================================================
        // 阶段 2：召唤 Mesh Optimizer 神库进行残暴切分！
        // ====================================================================
        const size_t max_vertices = 64;
        const size_t max_triangles = 124; // 硬件推荐 126，124 对齐更好
        const float cone_weight = 0.0f;   // 暂不考虑法线圆锥剔除优化

        // 询问神库：最多会切出多少个 Meshlet？
        size_t max_meshlets = meshopt_buildMeshletsBound(globalIndices.size(), max_vertices, max_triangles);

        // 准备临时的接收容器
        std::vector<meshopt_Meshlet> localMeshlets(max_meshlets);
        std::vector<unsigned int> meshletVertices(max_meshlets * max_vertices);
        std::vector<unsigned char> meshletTriangles(max_meshlets * max_triangles * 3);

        // ★ 执行切分！
        size_t meshlet_count = meshopt_buildMeshlets(
            localMeshlets.data(), meshletVertices.data(), meshletTriangles.data(),
            globalIndices.data(), globalIndices.size(), 
            &globalVertices[0].pos[0], globalVertices.size(), sizeof(MeshletVertex), 
            max_vertices, max_triangles, cone_weight);

        localMeshlets.resize(meshlet_count); // 裁剪掉多余的空间

        // ====================================================================
        // 阶段 3：重构数据，打造对 GPU 极度友好的内存布局
        // ====================================================================
        m_vertices.clear();
        m_microIndices.clear();
        m_meshlets.clear();

        for (const auto& m : localMeshlets) {
            MeshletDescriptor descriptor{};
            descriptor.vertexOffset = static_cast<uint32_t>(m_vertices.size());
            descriptor.triangleOffset = static_cast<uint32_t>(m_microIndices.size());
            descriptor.vertexCount = m.vertex_count;
            descriptor.triangleCount = m.triangle_count;
            m_meshlets.push_back(descriptor);

            // 1. 抽取这个 Meshlet 独占的物理顶点
            for (unsigned int i = 0; i < m.vertex_count; ++i) {
                m_vertices.push_back(globalVertices[meshletVertices[m.vertex_offset + i]]);
            }

            // 2. 抽取这个 Meshlet 的微型索引 (0~63 之间的数字)
            for (unsigned int i = 0; i < m.triangle_count * 3; ++i) {
                m_microIndices.push_back(static_cast<uint32_t>(meshletTriangles[m.triangle_offset + i]));
            }
        }

        std::cout << "[MeshletBuilder] SUCCESS! Model sliced into " << m_meshlets.size() << " highly optimized Meshlets!" << std::endl;
        std::cout << "[MeshletBuilder] Final buffers: Vertices(" << m_vertices.size() << "), MicroIndices(" << m_microIndices.size() << ")" << std::endl;

        return true;
    }

} // namespace Lizeral