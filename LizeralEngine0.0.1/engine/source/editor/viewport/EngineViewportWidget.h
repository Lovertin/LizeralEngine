#pragma once
#include "editor/EditorHeader.h"

class EngineViewportWidget : public QOpenGLWidget {
public:

    std::function<void()> onInitGL;
    
    EngineViewportWidget(Lizeral::Registry* registry, Lizeral::RenderingSystem* renderSys, QWidget* parent = nullptr) 
        : QOpenGLWidget(parent), m_Registry(registry), m_RenderSystem(renderSys) {


        // 允许此 Widget 接收键盘和鼠标事件
        setFocusPolicy(Qt::StrongFocus);
        setMouseTracking(true); 
    }

    void SetPhysicsSystem(Lizeral::PhysicsSystem* physSys) { m_physicsSystem = physSys; }

protected:
    // 1. OpenGL 上下文准备好时执行一次
    void initializeGL() override;

    // 2. 窗口大小改变时
    void resizeGL(int w, int h) override ;

    // 3. 每帧渲染逻辑
    void paintGL() override;

    // ---------------------------------------------------------
    // 【预留接口】：将 Qt 的鼠标键盘事件转发给引擎 Input 系统
    // ---------------------------------------------------------
    void mousePressEvent(QMouseEvent *event) override ;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override ;
    void keyPressEvent(QKeyEvent *event) override ;
    void keyReleaseEvent(QKeyEvent *event) override ;

private:
    Lizeral::Registry* m_Registry { nullptr };
    Lizeral::RenderingSystem* m_RenderSystem { nullptr };
    Lizeral::PhysicsSystem* m_physicsSystem {nullptr};

    bool m_isRoaming { false };
    QPoint m_roamStartGlobalPos; 
    QPointF m_virtualMousePos;
};
