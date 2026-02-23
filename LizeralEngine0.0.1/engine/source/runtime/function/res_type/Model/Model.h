#pragma once
#include <string>
#include <vector>
#include <memory>

#include "runtime/resource/resource.h"
#include "runtime/function/res_type/model/Mesh.h" // 请确保路径正确

// 引入 Assimp 核心头文件
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

namespace Lizeral {

    class Model : public Resource {
    public:
        Model() = default;
        virtual ~Model() = default;

        // 核心：实现 ResourceManager 的接口
        bool LoadFromFile(const std::string& path) override;

        // 渲染整个模型（遍历绘制所有的子网格）
        void Draw() const;

        // 获取内部所有的子网格数据
        std::vector<Mesh>& GetMeshes() { return m_Meshes; }
        const std::vector<Mesh>& GetMeshes() const { return m_Meshes; }

    private:
        std::vector<Mesh> m_Meshes;
        std::string m_Directory; // 保存模型所在的文件夹路径，为以后自动加载贴图做准备

        // 递归处理节点树
        void processNode(aiNode* node, const aiScene* scene);
        
        // 将 Assimp 的网格数据转化为 Lizeral 的网格数据
        Mesh processMesh(aiMesh* mesh, const aiScene* scene);
    };
}