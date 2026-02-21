#pragma once
#include "runtime/core/math/vector2.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/meta/reflection/reflection.h"
#include  "runtime/resource/resource.h"
#include <vector>
#include <memory>
#include<glad/glad.h>


namespace Lizeral{
    struct Vertex{
        Vector3 Position;
        Vector3 Normal;
        Vector2 TexCoords;

        // 比较操作符，用于可能的优化
        bool operator==(const Vertex& other) const {
            return Position == other.Position && 
                   Normal == other.Normal && 
                   TexCoords == other.TexCoords;
        }
    };
    // 可以添加颜色、切线等属性
    // Vector3 Color;
    // Vector3 Tangent;
    // Vector3 Bitangent;

    REFLECTION_TYPE(Mesh)
    CLASS(Mesh:public Resource,WhiteListFields) {
        REFLECTION_BODY(Mesh)
    public:
        Mesh() = default;
        Mesh(const std::vector<Vertex>& vertices, const std::vector<unsigned int>& indices);

        virtual ~Mesh();

        // 禁止拷贝，允许移动
        Mesh(const Mesh&) = delete;
        Mesh& operator=(const Mesh&) = delete;
        
        Mesh(Mesh&& other) noexcept;
        Mesh& operator=(Mesh&& other) noexcept;

        //载入Resource的虚函数
        virtual bool LoadFromFile(const std::string& path) override;

        //核心功能
        void Draw() const;  // const因为绘制不改变mesh状态
        void DrawInstanced(int instanceCount) const;
         // 获取器
        const std::vector<Vertex>& GetVertices() const { return m_Vertices; }
        const std::vector<unsigned int>& GetIndices() const { return m_Indices; }
        unsigned int GetVAO() const { return m_VAO; }
        size_t GetVertexCount() const { return m_Vertices.size(); }
        size_t GetIndexCount() const { return m_Indices.size(); }
        
        // 更新功能（谨慎使用，会重新上传数据到GPU）
        void UpdateVertices(const std::vector<Vertex>& vertices);
        void UpdateIndices(const std::vector<unsigned int>& indices);
        
        // 边界计算
        void CalculateBoundingBox(Vector3& min, Vector3& max) const;
        
        // 调试信息
        void PrintInfo() const;


    private:
        // 数据（CPU端）
        std::vector<Vertex> m_Vertices;
        std::vector<unsigned int> m_Indices;
        
        // OpenGL对象（GPU端）
        unsigned int m_VAO = 0;
        unsigned int m_VBO = 0;
        unsigned int m_EBO = 0;
        
        // 边界
        Vector3 m_BoundingBoxMin;
        Vector3 m_BoundingBoxMax;
        bool m_BoundingBoxValid = false;
        
    private:
        void SetupMesh();
        void UploadData(bool uploadVertices = true, bool uploadIndices = true);
        void UpdateBoundingBox();
        void ReleaseResources(); // 清理OpenGL资源

    };
}