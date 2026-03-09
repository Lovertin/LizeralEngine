#include "EngineViewportWidget.h"

void EngineViewportWidget::initializeGL(){
    // 【完美修复】：使用无捕获的 Lambda 包装 Qt 的 getProcAddress 成员函数
    if (!gladLoadGLLoader((GLADloadproc)[](const char* name) -> void* {
        QOpenGLContext* context = QOpenGLContext::currentContext();
        if (context) {
            // Qt 返回的是 QFunctionPointer，我们把它强转为 void* 喂给 GLAD
            return reinterpret_cast<void*>(context->getProcAddress(name));
        }
        return nullptr;
    })) {
        std::cout << "[Error] Failed to initialize GLAD with Qt Context!" << std::endl;
        return;
    }

    // GLAD 就绪后，初始化渲染系统
    if (m_RenderSystem) {
        m_RenderSystem->Initialize();
    }

    glEnable(GL_DEPTH_TEST);

    if (onInitGL) {
        onInitGL();
    }
}


void EngineViewportWidget::paintGL() {
    GLint qtCurrentFBO;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &qtCurrentFBO);

    float dpiScale = this->devicePixelRatio();

    QWidget* mainWindow = this->topLevelWidget(); 
    int fullW = static_cast<int>(mainWindow->width() * dpiScale);
    int fullH = static_cast<int>(mainWindow->height() * dpiScale);
    if (fullH == 0) fullH = 1;

    // 强制相机使用主窗口的物理比例
    float aspect = static_cast<float>(fullW) / static_cast<float>(fullH);
    if (m_Registry) {
        auto view = m_Registry->view<Lizeral::CameraComponent>();
        for (auto entity : view) {
            auto& cam = view.get<Lizeral::CameraComponent>(entity);
            if (cam.isMain()) { 
                cam.setAspect(aspect);
                break;
            }
        }
    }

    QPoint topLeftPos = this->mapTo(mainWindow, QPoint(0, 0));
    
    // 把 Qt 的逻辑坐标转换成 OpenGL 需要的物理坐标！
    int physicalOffsetX = static_cast<int>(topLeftPos.x() * dpiScale);
    
    int physicalWidgetH = static_cast<int>(this->height() * dpiScale);
    int physicalOffsetY = fullH - static_cast<int>(topLeftPos.y() * dpiScale + physicalWidgetH);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (m_RenderSystem && m_Registry) {
        m_RenderSystem->SetDefaultFBO(static_cast<unsigned int>(qtCurrentFBO));
        
        // 传入绝对精准的物理坐标抵消！
        m_RenderSystem->SetViewport(-physicalOffsetX, -physicalOffsetY, fullW, fullH);

        m_RenderSystem->Tick(*m_Registry, 0.016f); 
    }

    if (m_physicsSystem && m_Registry) {
        // 从 ECS 里抓取主摄像机的数据
        Lizeral::Matrix4x4 viewMatrix = Lizeral::Matrix4x4::IDENTITY;
        Lizeral::Matrix4x4 projMatrix = Lizeral::Matrix4x4::IDENTITY;

        auto cameraView = m_Registry->view<Lizeral::CameraComponent>();
        for (auto entity : cameraView) {
            auto& cam = cameraView.get<Lizeral::CameraComponent>(entity);
            if (cam.isMain()) {
                viewMatrix = cam.getViewMatrix();
                projMatrix = cam.getProjectionMatrix();
                break;
            }
        }

        glUseProgram(0);

        m_physicsSystem->DebugDrawWorld();
        
        // 2. 获取真正的 Drawer 实例，并调用现代的批量绘制函数！
        auto* debugDrawer = static_cast<Lizeral::PhysicsDebugDrawer*>(m_physicsSystem->getScene()->getWorld()->getDebugDrawer());
        if (debugDrawer) {
            debugDrawer->FlushLines(viewMatrix, projMatrix);
        }
    }
}

void EngineViewportWidget::resizeGL(int w, int h)
{
    // glViewport(0, 0, w, h);
}

void EngineViewportWidget::mousePressEvent(QMouseEvent *event){
    if (event->button() == Qt::RightButton) {
        m_isRoaming = true;
        
        // 1. 记录右键按下时的【系统全局坐标】（锚点）
        m_roamStartGlobalPos = QCursor::pos(); 
        
        // 2. 将当前的局部坐标作为虚拟坐标的起点
        m_virtualMousePos = event->pos(); 

        Lizeral::Input::GetInstance().SetMouseButtonDown(Lizeral::MouseButton::Right, true);
        
        // 3. 隐藏光标
        setCursor(Qt::BlankCursor);
    }
}

void EngineViewportWidget::mouseReleaseEvent(QMouseEvent *event){
    if (event->button() == Qt::RightButton) {
        m_isRoaming = false;

        Lizeral::Input::GetInstance().SetMouseButtonDown(Lizeral::MouseButton::Right, false);
        setCursor(Qt::ArrowCursor); 
        QCursor::setPos(m_roamStartGlobalPos);
    }
}

void EngineViewportWidget::mouseMoveEvent(QMouseEvent *event){
    if (m_isRoaming) {

        QPoint currentGlobalPos = QCursor::pos();

        if (currentGlobalPos == m_roamStartGlobalPos) {
            return;
        }

        int dx = currentGlobalPos.x() - m_roamStartGlobalPos.x();
        int dy = currentGlobalPos.y() - m_roamStartGlobalPos.y();

        m_virtualMousePos.setX(m_virtualMousePos.x() + dx);
        m_virtualMousePos.setY(m_virtualMousePos.y() + dy);

        Lizeral::Input::GetInstance().SetMousePosition(m_virtualMousePos.x(), m_virtualMousePos.y());

        QCursor::setPos(m_roamStartGlobalPos);
        
    } else {
        // 普通模式下，正常传递坐标给引擎即可（用于后续的 UI 点击、射线拾取等）
        Lizeral::Input::GetInstance().SetMousePosition(event->pos().x(), event->pos().y());
    }
}

void EngineViewportWidget::keyPressEvent(QKeyEvent *event){
    if (event->key() == Qt::Key_W) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::W, true);
    if (event->key() == Qt::Key_S) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::S, true);
    if (event->key() == Qt::Key_A) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::A, true);
    if (event->key() == Qt::Key_D) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::D, true);
    if (event->key() == Qt::Key_Q) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::Q, true);
    if (event->key() == Qt::Key_E) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::E, true);
}

void EngineViewportWidget::keyReleaseEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_W) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::W, false);
    if (event->key() == Qt::Key_S) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::S, false);
    if (event->key() == Qt::Key_A) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::A, false);
    if (event->key() == Qt::Key_D) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::D, false);
    if (event->key() == Qt::Key_Q) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::Q, false);
    if (event->key() == Qt::Key_E) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::E, false);
}