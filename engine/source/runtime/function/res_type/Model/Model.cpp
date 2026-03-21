#include "Model.h"
#include "runtime/function/res_type/Material/PBRMaterial.h"
#include "runtime/resource/resourceManager/resourceManager.h"
#include "runtime/function/res_type/texture/Texture2D.h"
#include <iostream>

namespace Lizeral {

    bool Model::LoadFromFile(const std::string& path) {
        Assimp::Importer importer;
        std::cout << "[Assimp] Ready to call ReadFile for: " << path << std::endl;

        const aiScene* scene = importer.ReadFile(path, 
            aiProcess_Triangulate | 
            aiProcess_GenSmoothNormals | 
            aiProcess_FlipUVs | 
            aiProcess_CalcTangentSpace
            // aiProcess_PreTransformVertices
        );

        std::cout << "[Assimp] ReadFile completed!" << std::endl;

        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "[Assimp] ERROR:: " << importer.GetErrorString() << std::endl;
            return false;
        }

        //  "asset/model/shuttle/"
        m_Directory = path.substr(0, path.find_last_of('/') != std::string::npos ? path.find_last_of('/') : path.find_last_of('\\'));

        processNode(scene->mRootNode, scene);

        std::cout << "[Model] Successfully loaded model: " << path << " (Sub-meshes: " << m_Meshes.size() << ")" << std::endl;
        return true;
    }

    void Model::Draw() const {
        for (unsigned int i = 0; i < m_Meshes.size(); i++) {
            m_Meshes[i].Draw();
        }
    }

    void Model::processNode(aiNode* node, const aiScene* scene) {
        for (unsigned int i = 0; i < node->mNumMeshes; i++) {
            // indices node->mMeshes, data scene->mMeshes
            aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
            m_Meshes.push_back(processMesh(mesh, scene));
        }

        for (unsigned int i = 0; i < node->mNumChildren; i++) {
            processNode(node->mChildren[i], scene);
        }
    }

    Mesh Model::processMesh(aiMesh* mesh, const aiScene* scene) {
        std::vector<Vertex> vertices;
        std::vector<unsigned int> indices;

        for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
            Vertex vertex;

            vertex.Position = Vector3(mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z);

            if (mesh->HasNormals()) {
                vertex.Normal = Vector3(mesh->mNormals[i].x, mesh->mNormals[i].y, mesh->mNormals[i].z);
            }

            if (mesh->mTextureCoords[0]) {
                vertex.TexCoords = Vector2(mesh->mTextureCoords[0][i].x, mesh->mTextureCoords[0][i].y);

                if (mesh->HasTangentsAndBitangents()) {
                    vertex.Tangent = Vector3(mesh->mTangents[i].x, mesh->mTangents[i].y, mesh->mTangents[i].z);
                    vertex.Bitangent = Vector3(mesh->mBitangents[i].x, mesh->mBitangents[i].y, mesh->mBitangents[i].z);
                }
            } else {
                vertex.TexCoords = Vector2(0.0f, 0.0f);
            }

            vertices.push_back(vertex);
        }

        for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
            aiFace face = mesh->mFaces[i];
            for (unsigned int j = 0; j < face.mNumIndices; j++) {
                indices.push_back(face.mIndices[j]);
            }
        }

        auto pbrMat = std::make_shared<PBRMaterial>();

        if (mesh->mMaterialIndex >= 0) {
            aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
            
            // Base Color / Diffuse Color
            aiColor4D diffuseColor;
            // Assimp put glTF baseColorFactor to AI_MATKEY_COLOR_DIFFUSE
            if (AI_SUCCESS == aiGetMaterialColor(material, AI_MATKEY_COLOR_DIFFUSE, &diffuseColor)) {
                pbrMat->m_Albedo = Vector3(diffuseColor.r, diffuseColor.g, diffuseColor.b);
            }

            // read Albedo Map
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