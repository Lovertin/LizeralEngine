#include "Model.h"
#include "runtime/function/res_type/Material/PBRMaterial.h"
#include "runtime/resource/resourceManager/resourceManager.h"
#include "runtime/function/res_type/texture/Texture2D.h"
#include <iostream>

namespace Lizeral {

    bool Model::LoadFromFile(const std::string& path) {
        Assimp::Importer importer;
        std::cout << "[Assimp] Ready to call ReadFile for: " << path << std::endl;

        // 【超级魔法宏】：
        // aiProcess_Triangulate: 如果模型有四边形面，全部切成三角形（完美解决你的破面问题！）
        // aiProcess_GenSmoothNormals: 如果模型没有法线，自动平滑计算法线
        // aiProcess_FlipUVs: OpenGL 的 UV 原点在左下角，自动帮你翻转过来
        // aiProcess_CalcTangentSpace: 自动计算切线和副切线（为以后的高品质法线贴图做准备）
        const aiScene* scene = importer.ReadFile(path, 
            aiProcess_Triangulate | 
            aiProcess_GenSmoothNormals | 
            aiProcess_FlipUVs | 
            aiProcess_CalcTangentSpace
            // aiProcess_PreTransformVertices
        );

        std::cout << "[Assimp] ReadFile completed!" << std::endl;

        // 错误检查
        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "[Assimp] ERROR:: " << importer.GetErrorString() << std::endl;
            return false;
        }

        // 保存文件目录路径，比如 "asset/model/shuttle/"
        m_Directory = path.substr(0, path.find_last_of('/') != std::string::npos ? path.find_last_of('/') : path.find_last_of('\\'));

        // 开始递归处理节点树
        processNode(scene->mRootNode, scene);

        std::cout << "[Model] Successfully loaded model: " << path << " (Sub-meshes: " << m_Meshes.size() << ")" << std::endl;
        return true;
    }

    void Model::Draw() const {
        // 遍历渲染所有的子网格
        for (unsigned int i = 0; i < m_Meshes.size(); i++) {
            m_Meshes[i].Draw();
        }
    }

    void Model::processNode(aiNode* node, const aiScene* scene) {
        // 1. 处理当前节点挂载的所有网格
        for (unsigned int i = 0; i < node->mNumMeshes; i++) {
            // node->mMeshes 保存的是索引，真正的数据在 scene->mMeshes 里
            aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
            m_Meshes.push_back(processMesh(mesh, scene));
        }

        // 2. 递归处理所有的子节点
        for (unsigned int i = 0; i < node->mNumChildren; i++) {
            processNode(node->mChildren[i], scene);
        }
    }

    Mesh Model::processMesh(aiMesh* mesh, const aiScene* scene) {
        std::vector<Vertex> vertices;
        std::vector<unsigned int> indices;

        // ==========================================
        // 1. 提取顶点数据 (Vertices)
        // ==========================================
        for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
            Vertex vertex;
            
            // 位置
            vertex.Position = Vector3(mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z);
            
            // 法线
            if (mesh->HasNormals()) {
                vertex.Normal = Vector3(mesh->mNormals[i].x, mesh->mNormals[i].y, mesh->mNormals[i].z);
            }

            // 纹理坐标 (UV) 和 切线空间
            // Assimp 允许一个顶点有 8 套 UV，我们目前只关心第一套 (索引 0)
            if (mesh->mTextureCoords[0]) {
                vertex.TexCoords = Vector2(mesh->mTextureCoords[0][i].x, mesh->mTextureCoords[0][i].y);
                
                // 切线
                if (mesh->HasTangentsAndBitangents()) {
                    vertex.Tangent = Vector3(mesh->mTangents[i].x, mesh->mTangents[i].y, mesh->mTangents[i].z);
                    vertex.Bitangent = Vector3(mesh->mBitangents[i].x, mesh->mBitangents[i].y, mesh->mBitangents[i].z);
                }
            } else {
                vertex.TexCoords = Vector2(0.0f, 0.0f);
            }

            vertices.push_back(vertex);
        }

        // ==========================================
        // 2. 提取索引数据 (Indices)
        // ==========================================
        for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
            aiFace face = mesh->mFaces[i];
            // 之前我们在 ReadFile 时传了 aiProcess_Triangulate，所以这里的 mNumIndices 必然是 3
            for (unsigned int j = 0; j < face.mNumIndices; j++) {
                indices.push_back(face.mIndices[j]);
            }
        }

        // ==========================================
        // 3. 材质处理 (暂时传入 nullptr，解耦渲染)
        // ==========================================
        // 真正的商业引擎会在这里读取 Assimp 的材质，并生成对应的 PBRMaterial。
        // 但为了保证本次重构的安全着陆，我们先不自动加载贴图，
        // 我们返回 nullptr，允许你在外部统一绑定 Shader。
        auto pbrMat = std::make_shared<PBRMaterial>();

        if (mesh->mMaterialIndex >= 0) {
            aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
            
            // --- A. 尝试读取纯色 (Base Color / Diffuse Color) ---
            aiColor4D diffuseColor;
            // Assimp 会将 glTF 的 baseColorFactor 映射到 AI_MATKEY_COLOR_DIFFUSE
            if (AI_SUCCESS == aiGetMaterialColor(material, AI_MATKEY_COLOR_DIFFUSE, &diffuseColor)) {
                // 将读取到的颜色赋值给材质的 u_Albedo (假设你的 PBRMaterial 叫 m_Albedo 或 m_Color)
                // 注意：glTF 里的颜色通常已经是线性的了，直接赋值即可
                pbrMat->m_Albedo = Vector3(diffuseColor.r, diffuseColor.g, diffuseColor.b);
            }

            // --- B. 尝试读取贴图 (Albedo Map) 会覆盖上面的纯色 ---
            aiString str;
            if (material->GetTextureCount(aiTextureType_DIFFUSE) > 0) {
                material->GetTexture(aiTextureType_DIFFUSE, 0, &str);
                
                const aiTexture* embeddedTex = scene->GetEmbeddedTexture(str.C_Str());
                if (embeddedTex) {
                    auto tex = std::make_shared<Texture2D>();
                    if (embeddedTex->mHeight == 0) {
                        tex->LoadFromMemory(reinterpret_cast<unsigned char*>(embeddedTex->pcData), embeddedTex->mWidth);
                    }
                    pbrMat->m_AlbedoMap = tex;
                } else {
                    std::string texPath = m_Directory + "\\" + str.C_Str();
                    pbrMat->m_AlbedoMap = ResourceManager::GetInstance().Load<Texture2D>(texPath);
                }
            }
        }

        return Mesh(vertices, indices, pbrMat);
    }
}