#include "PhysicsDebugDrawer.h"
#include <glad/glad.h> // 确保引入了你的 OpenGL 头文件
#include <iostream>

namespace Lizeral {

    PhysicsDebugDrawer::PhysicsDebugDrawer() {
        // 默认开启绘制线框 (Wireframe)
        m_debugMode = btIDebugDraw::DBG_DrawWireframe; 
    }

    // 【核心魔法】：Bullet 每画一条边，都会调用这个函数！
    void PhysicsDebugDrawer::drawLine(const btVector3& from, const btVector3& to, const btVector3& color) {
        // 保存 OpenGL 状态
        glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT);
        
        glDisable(GL_DEPTH_TEST); // 关闭深度测试，让你能透过墙壁和模型看到物理线框
        glDisable(GL_LIGHTING);   // 关闭光照，纯色显示

        glDisable(GL_TEXTURE_2D); 
        glDisable(GL_BLEND);

        glLineWidth(2.0f);        // 线条加粗

        // 使用传统的立即模式画线
        glBegin(GL_LINES);
        glColor3f(color.getX(), color.getY(), color.getZ()); // Bullet 会自动分配颜色！静态物体是绿的，动态是白的，休眠是红的！
        glVertex3f(from.getX(), from.getY(), from.getZ());
        glVertex3f(to.getX(), to.getY(), to.getZ());
        glEnd();

        glPopAttrib(); // 恢复状态
    }

    // 其他几个纯虚函数我们可以先空着，但必须实现
    void PhysicsDebugDrawer::drawContactPoint(const btVector3& PointOnB, const btVector3& normalOnB, btScalar distance, int lifeTime, const btVector3& color) {}
    void PhysicsDebugDrawer::reportErrorWarning(const char* warningString) { std::cerr << "[Physics Debug] " << warningString << std::endl; }
    void PhysicsDebugDrawer::draw3dText(const btVector3& location, const char* textString) {}
    
    void PhysicsDebugDrawer::setDebugMode(int debugMode) { m_debugMode = debugMode; }
    int PhysicsDebugDrawer::getDebugMode() const { return m_debugMode; }

} // namespace Lizeral