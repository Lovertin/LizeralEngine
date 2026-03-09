#include "PhysicsDebugDrawer.h"
#include <glad/glad.h>
#include <iostream>

namespace Lizeral {

    PhysicsDebugDrawer::PhysicsDebugDrawer() {
        m_debugMode = btIDebugDraw::DBG_DrawWireframe; 
    }

    PhysicsDebugDrawer::~PhysicsDebugDrawer() {
        if (m_VAO != 0) glDeleteVertexArrays(1, &m_VAO);
        if (m_VBO != 0) glDeleteBuffers(1, &m_VBO);
    }

    // 【极速操作】：只把顶点和颜色塞进数组，绝对不在这里调用 gl 命令！
    void PhysicsDebugDrawer::drawLine(const btVector3& from, const btVector3& to, const btVector3& color) {
        // 起点位置 + 颜色
        m_lineData.push_back(from.x()); m_lineData.push_back(from.y()); m_lineData.push_back(from.z());
        m_lineData.push_back(color.x()); m_lineData.push_back(color.y()); m_lineData.push_back(color.z());
        
        // 终点位置 + 颜色
        m_lineData.push_back(to.x()); m_lineData.push_back(to.y()); m_lineData.push_back(to.z());
        m_lineData.push_back(color.x()); m_lineData.push_back(color.y()); m_lineData.push_back(color.z());
    }

    // 【大招释放】：一帧调用一次，把刚才攒的几千个点，一口气画出来！
    void PhysicsDebugDrawer::FlushLines(const Matrix4x4& viewMatrix, const Matrix4x4& projMatrix) {
        if (m_lineData.empty()) return;

        // 1. 懒加载：初始化专属的 Line Shader 和 VAO/VBO
        if (m_VAO == 0) {
            // 注意：请将这里的路径替换为你电脑里的实际路径！
            m_lineShader = std::make_shared<Shader>(
                "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\line.vert",
                "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\line.frag"
            );

            glGenVertexArrays(1, &m_VAO);
            glGenBuffers(1, &m_VBO);

            glBindVertexArray(m_VAO);
            glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
            // 属性0：位置 (vec3)
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
            // 属性1：颜色 (vec3)
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
            glBindVertexArray(0);
        }

        // 2. 绑定 Shader 并上传矩阵
        m_lineShader->Bind();
        m_lineShader->SetUniformMat4f("view", viewMatrix);
        m_lineShader->SetUniformMat4f("projection", projMatrix);

        // 3. 将本帧收集的所有线段数据上传到显存 (使用 GL_DYNAMIC_DRAW 因为每帧都在变)
        glBindVertexArray(m_VAO);
        glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
        glBufferData(GL_ARRAY_BUFFER, m_lineData.size() * sizeof(float), m_lineData.data(), GL_DYNAMIC_DRAW);

        // 4. 开始绘制！
        glDisable(GL_DEPTH_TEST); // 关闭深度测试，实现透视查看
        glLineWidth(2.0f);        // 线条加粗
        
        // 核心：画线！顶点数量 = 数组总长度 / 6 (因为一个点占 6 个 float)
        glDrawArrays(GL_LINES, 0, m_lineData.size() / 6); 

        glEnable(GL_DEPTH_TEST);  // 恢复状态
        glBindVertexArray(0);
        m_lineShader->Unbind();

        // 5. 【极其重要】：清空数组，迎接下一帧！
        m_lineData.clear();
    }

    void PhysicsDebugDrawer::drawContactPoint(const btVector3& PointOnB, const btVector3& normalOnB, btScalar distance, int lifeTime, const btVector3& color) {}
    void PhysicsDebugDrawer::reportErrorWarning(const char* warningString) { std::cerr << "[Physics Debug] " << warningString << std::endl; }
    void PhysicsDebugDrawer::draw3dText(const btVector3& location, const char* textString) {}
    void PhysicsDebugDrawer::setDebugMode(int debugMode) { m_debugMode = debugMode; }
    int PhysicsDebugDrawer::getDebugMode() const { return m_debugMode; }

} // namespace Lizeral