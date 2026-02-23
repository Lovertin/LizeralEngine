#pragma once

#include <LinearMath/btIDebugDraw.h>
// 引入你的数学库或者通用头文件
#include "runtime/core/math/vector3.h"

namespace Lizeral {

    class PhysicsDebugDrawer : public btIDebugDraw {
    private:
        int m_debugMode;

    public:
        PhysicsDebugDrawer();
        virtual ~PhysicsDebugDrawer() = default;

        // --- 必须重写的核心绘制函数 ---
        virtual void drawLine(const btVector3& from, const btVector3& to, const btVector3& color) override;
        virtual void drawContactPoint(const btVector3& PointOnB, const btVector3& normalOnB, btScalar distance, int lifeTime, const btVector3& color) override;
        virtual void reportErrorWarning(const char* warningString) override;
        virtual void draw3dText(const btVector3& location, const char* textString) override;
        virtual void setDebugMode(int debugMode) override;
        virtual int getDebugMode() const override;
    };

} // namespace Lizeral