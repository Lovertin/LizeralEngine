#pragma once

#include <vector>
#include <memory>
#include <glad/glad.h>
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/vector2.h"
#include "runtime/function/res_type/Material/Material.h" 
#include "runtime/core/meta/reflection/reflection.h"

namespace Lizeral {

    struct Vertex {
        Vector3 Position;
        Vector3 Normal;
        Vector2 TexCoords;
        Vector3 Tangent;   
        Vector3 Bitangent; 
    };
    REFLECTION_TYPE(Mesh)
    CLASS(Mesh,WhiteListFields) {
        REFLECTION_BODY(Mesh)
    public:
        std::vector<Vertex>       m_Vertices;
        std::vector<unsigned int> m_Indices;

        std::shared_ptr<Material> m_Material;

        Mesh() = default;

        Mesh(std::vector<Vertex> vertices, 
             std::vector<unsigned int> indices, 
             std::shared_ptr<Material> material);

        void Draw() const;

    private:
        unsigned int VAO, VBO, EBO;
        void setupMesh(); 
    };
}