#pragma once
#include <LinearMath/btIDebugDraw.h>
#include <vector>
#include <memory>
#include "runtime/function/res_type/shader/shader.h" // 引入你的 Shader 类
#include "runtime/core/math/matrix4.h"             // 引入你的矩阵类

namespace Lizeral {

    class PhysicsDebugDrawer : public btIDebugDraw {
    public:
        PhysicsDebugDrawer();
        virtual ~PhysicsDebugDrawer();

        virtual void drawLine(const btVector3& from, const btVector3& to, const btVector3& color) override;
        virtual void drawContactPoint(const btVector3& PointOnB, const btVector3& normalOnB, btScalar distance, int lifeTime, const btVector3& color) override;
        virtual void reportErrorWarning(const char* warningString) override;
        virtual void draw3dText(const btVector3& location, const char* textString) override;
        virtual void setDebugMode(int debugMode) override;
        virtual int getDebugMode() const override;

        void FlushLines(const Matrix4x4& viewMatrix, const Matrix4x4& projMatrix);

    private:
        int m_debugMode;
        
        // 数据缓存：按照 [X, Y, Z, R, G, B] 的顺序存入
        std::vector<float> m_lineData;

        // 现代 OpenGL 句柄
        unsigned int m_VAO = 0;
        unsigned int m_VBO = 0;
        std::shared_ptr<Shader> m_lineShader = nullptr;
    };

} // namespace Lizeral