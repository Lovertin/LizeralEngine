#pragma once

#include <vector>
#include <memory>
#include <glad/glad.h>
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/vector2.h"
#include "runtime/function/res_type/Material/Material.h" // 引入材质基类
#include "runtime/core/meta/reflection/reflection.h"

namespace Lizeral {

    // 工业级顶点结构体（内存连续，完美契合 OpenGL）
    struct Vertex {
        Vector3 Position;
        Vector3 Normal;
        Vector2 TexCoords;
        Vector3 Tangent;   // 【新增】：切线！用于法线贴图 (Normal Mapping) 极其重要
        Vector3 Bitangent; // 【新增】：副切线
    };
    REFLECTION_TYPE(Mesh)
    CLASS(Mesh,WhiteListFields) {
        REFLECTION_BODY(Mesh)
    public:
        // 网格数据
        std::vector<Vertex>       m_Vertices;
        std::vector<unsigned int> m_Indices;
        
        // 【关键】：每个子网格持有自己特有的材质实例
        std::shared_ptr<Material> m_Material;

        Mesh() = default;

        // 构造函数：接收数据并立刻传给显卡
        Mesh(std::vector<Vertex> vertices, 
             std::vector<unsigned int> indices, 
             std::shared_ptr<Material> material);

        // 渲染这个子网格
        void Draw() const;

    private:
        // 渲染数据 ID
        unsigned int VAO, VBO, EBO;
        void setupMesh(); // 内部初始化 OpenGL 缓冲
    };
}