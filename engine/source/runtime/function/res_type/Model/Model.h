#pragma once
#include <string>
#include <vector>
#include <memory>

#include "runtime/resource/resource.h"
#include "runtime/function/res_type/model/Mesh.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

namespace Lizeral {

    class Model : public Resource {
    public:
        Model() = default;
        virtual ~Model() = default;

        bool LoadFromFile(const std::string& path) override;

        void Draw() const;

        std::vector<Mesh>& GetMeshes() { return m_Meshes; }
        const std::vector<Mesh>& GetMeshes() const { return m_Meshes; }

    private:
        std::vector<Mesh> m_Meshes;
        std::string m_Directory; 

        void processNode(aiNode* node, const aiScene* scene);
        
        Mesh processMesh(aiMesh* mesh, const aiScene* scene);
    };
}