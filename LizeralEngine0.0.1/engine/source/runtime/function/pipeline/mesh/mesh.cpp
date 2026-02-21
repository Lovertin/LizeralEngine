// runtime/core/mesh.cpp
#include "mesh.h"
#include "runtime/function/pipeline/utils/ImportModel.h"
#include <iostream>
#include <limits>

namespace Lizeral {

    Mesh::Mesh(const std::vector<Vertex>& vertices, const std::vector<unsigned int>& indices)
        : m_Vertices(vertices), m_Indices(indices) {
        
        if (m_Vertices.empty()) {
            std::cerr << "Warning: Creating mesh with no vertices!" << std::endl;
            return;
        }
        
        if (m_Indices.empty()) {
            std::cerr << "Warning: Creating mesh with no indices!" << std::endl;
        }
        
        SetupMesh();
        UpdateBoundingBox();
    }

    Mesh::~Mesh() {
        ReleaseResources();
    }

    Mesh::Mesh(Mesh&& other) noexcept
        : m_Vertices(std::move(other.m_Vertices))
        , m_Indices(std::move(other.m_Indices))
        , m_VAO(other.m_VAO)
        , m_VBO(other.m_VBO)
        , m_EBO(other.m_EBO)
        , m_BoundingBoxMin(other.m_BoundingBoxMin)
        , m_BoundingBoxMax(other.m_BoundingBoxMax)
        , m_BoundingBoxValid(other.m_BoundingBoxValid) {
        
        // 清空源对象的OpenGL句柄，防止其在析构时释放资源
        other.m_VAO = 0;
        other.m_VBO = 0;
        other.m_EBO = 0;
    }

    Mesh& Mesh::operator=(Mesh&& other) noexcept {
        if (this != &other) {
            // 释放当前资源
            ReleaseResources();
            
            // 移动数据
            m_Vertices = std::move(other.m_Vertices);
            m_Indices = std::move(other.m_Indices);
            m_VAO = other.m_VAO;
            m_VBO = other.m_VBO;
            m_EBO = other.m_EBO;
            m_BoundingBoxMin = other.m_BoundingBoxMin;
            m_BoundingBoxMax = other.m_BoundingBoxMax;
            m_BoundingBoxValid = other.m_BoundingBoxValid;
            
            // 清空源对象
            other.m_VAO = 0;
            other.m_VBO = 0;
            other.m_EBO = 0;
        }
        return *this;
    }

    bool Mesh::LoadFromFile(const std::string& path) {
        ImportModel loader(path.c_str());
        m_Vertices = loader.GetMeshVertices();
        m_Indices.clear(); // .obj 铺平后没有索引

        ReleaseResources(); // 清理旧数据

        // 【关键】：把数据从 CPU 内存打包装入 GPU 显存！
        SetupMesh(); 
        
        UpdateBoundingBox();
        return true;
    }

    void Mesh::UploadData(bool uploadVertices, bool uploadIndices) {
        glBindVertexArray(m_VAO);
        
        if (uploadVertices && !m_Vertices.empty()) {
            glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
            glBufferData(GL_ARRAY_BUFFER, m_Vertices.size() * sizeof(Vertex), 
                        m_Vertices.data(), GL_STATIC_DRAW);
        }
        
        if (uploadIndices && !m_Indices.empty()) {
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_Indices.size() * sizeof(unsigned int), 
                        m_Indices.data(), GL_STATIC_DRAW);
        }
        
        glBindVertexArray(0);
    }

    void Mesh::SetupMesh() {
        if (m_Vertices.empty()) return;

        // 生成 VAO 和 VBO 的唯一 ID
        glGenVertexArrays(1, &m_VAO);
        glGenBuffers(1, &m_VBO);

        // 绑定 VAO (开始记录状态)
        glBindVertexArray(m_VAO);

        // 绑定 VBO 并把 m_Vertices 的数据送入 GPU 显存！
        glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
        glBufferData(GL_ARRAY_BUFFER, m_Vertices.size() * sizeof(Vertex), &m_Vertices[0], GL_STATIC_DRAW);

        // 告诉显卡，这些数据是什么格式的：
        // a. 顶点坐标 Position (对应 Shader 里的 layout (location = 0))
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        
        // b. 法线 Normal (对应 Shader 里的 layout (location = 1))
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, Normal));
        
        // c. 纹理坐标 TexCoords (对应 Shader 里的 layout (location = 2))
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, TexCoords));

        // 如果有索引，也配置 EBO
        if (!m_Indices.empty()) {
            glGenBuffers(1, &m_EBO);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_Indices.size() * sizeof(unsigned int), &m_Indices[0], GL_STATIC_DRAW);
        }

        // 解绑 VAO (记录结束)
        glBindVertexArray(0);
        
        std::cout << "[Mesh] Successfully uploaded " << m_Vertices.size() << " vertices to GPU! VAO ID: " << m_VAO << std::endl;
    }

    void Mesh::Draw() const {
        if (m_VAO == 0) {
            // 🚨 如果你看到了这行输出，说明元凶找到了！
            std::cout << "[ERROR] Mesh::Draw() aborted! VAO is 0. Data was never uploaded to GPU!" << std::endl;
            return;
        }

        // 1. 绑定 VAO
        glBindVertexArray(m_VAO);

        // 2. 智能绘制判断
        if (!m_Indices.empty()) {
            // 如果有索引，使用 EBO 绘制 (Indexed Drawing)
            glDrawElements(GL_TRIANGLES, m_Indices.size(), GL_UNSIGNED_INT, 0);
        } else {
            // 【关键修复】：如果没有索引，直接顺序绘制顶点！(Array Drawing)
            glDrawArrays(GL_TRIANGLES, 0, m_Vertices.size());
        }

        // 3. 解绑 VAO（好习惯，防止状态污染）
        glBindVertexArray(0);
    }

    void Mesh::DrawInstanced(int instanceCount) const {
        if (m_VAO == 0 || m_Vertices.empty() || instanceCount <= 0) {
            return;
        }
        
        glBindVertexArray(m_VAO);
        
        if (!m_Indices.empty()) {
            glDrawElementsInstanced(GL_TRIANGLES, static_cast<GLsizei>(m_Indices.size()), 
                                   GL_UNSIGNED_INT, 0, instanceCount);
        } else {
            glDrawArraysInstanced(GL_TRIANGLES, 0, static_cast<GLsizei>(m_Vertices.size()), 
                                 instanceCount);
        }
        
        glBindVertexArray(0);
    }

    void Mesh::UpdateVertices(const std::vector<Vertex>& vertices) {
        m_Vertices = vertices;
        UploadData(true, false);  // 只更新顶点数据
        UpdateBoundingBox();
    }

    void Mesh::UpdateIndices(const std::vector<unsigned int>& indices) {
        m_Indices = indices;
        UploadData(false, true);  // 只更新索引数据
    }

    void Mesh::UpdateBoundingBox() {
        if (m_Vertices.empty()) {
            m_BoundingBoxValid = false;
            return;
        }
        
        const float INF = std::numeric_limits<float>::max();
        m_BoundingBoxMin = Vector3(INF, INF, INF);
        m_BoundingBoxMax = Vector3(-INF, -INF, -INF);
        
        for (const auto& vertex : m_Vertices) {
            m_BoundingBoxMin.x = std::min(m_BoundingBoxMin.x, vertex.Position.x);
            m_BoundingBoxMin.y = std::min(m_BoundingBoxMin.y, vertex.Position.y);
            m_BoundingBoxMin.z = std::min(m_BoundingBoxMin.z, vertex.Position.z);
            
            m_BoundingBoxMax.x = std::max(m_BoundingBoxMax.x, vertex.Position.x);
            m_BoundingBoxMax.y = std::max(m_BoundingBoxMax.y, vertex.Position.y);
            m_BoundingBoxMax.z = std::max(m_BoundingBoxMax.z, vertex.Position.z);
        }
        
        m_BoundingBoxValid = true;
    }

    void Mesh::CalculateBoundingBox(Vector3& min, Vector3& max) const {
        if (m_BoundingBoxValid) {
            min = m_BoundingBoxMin;
            max = m_BoundingBoxMax;
        } else {
            min = max = Vector3(0, 0, 0);
        }
    }

    void Mesh::PrintInfo() const {
        std::cout << "Mesh Information:" << std::endl;
        std::cout << "  Vertices: " << m_Vertices.size() << std::endl;
        std::cout << "  Indices: " << m_Indices.size() << std::endl;
        std::cout << "  VAO: " << m_VAO << std::endl;
        std::cout << "  VBO: " << m_VBO << std::endl;
        std::cout << "  EBO: " << m_EBO << std::endl;
        
        if (m_BoundingBoxValid) {
            std::cout << "  Bounding Box:" << std::endl;
            std::cout << "    Min: (" << m_BoundingBoxMin.x << ", " 
                      << m_BoundingBoxMin.y << ", " << m_BoundingBoxMin.z << ")" << std::endl;
            std::cout << "    Max: (" << m_BoundingBoxMax.x << ", " 
                      << m_BoundingBoxMax.y << ", " << m_BoundingBoxMax.z << ")" << std::endl;
        }
    }

    void Mesh::ReleaseResources() {
        if (m_VAO != 0) {
            glDeleteVertexArrays(1, &m_VAO);
            m_VAO = 0;
        }
        if (m_VBO != 0) {
            glDeleteBuffers(1, &m_VBO);
            m_VBO = 0;
        }
        if (m_EBO != 0) {
            glDeleteBuffers(1, &m_EBO);
            m_EBO = 0;
        }
    }

} // namespace Lizeral