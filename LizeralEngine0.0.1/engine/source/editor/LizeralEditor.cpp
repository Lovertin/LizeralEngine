#include <glad/glad.h>
// 注意：在 Qt 中，我们尽量不使用 GLFW，因为窗口和输入由 Qt 接管了
// 但如果你的底层代码(比如 Input)还强依赖它，可以保留。
#include <GLFW/glfw3.h>

#include <QApplication>
#include <QMainWindow>
#include <QDockWidget>
#include <QScrollArea>
#include <QOpenGLWidget>
#include <QOpenGLContext> // 【新增】：用于获取 Qt 的 OpenGL 函数指针
#include <QVBoxLayout>
#include <QLabel>
#include <iostream>
#include <QTimer>
#include <QElapsedTimer>
#include <QMouseEvent>
#include <QKeyEvent>

// --- 引擎侧模块 ---
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/componentAll.h"

// --- 编辑器侧模块 ---
#include "editor/factory/DataTypeUIFactory/DataTypeFactory.h"
#include "editor/selection/EditorSelection.h"
#include "editor/panels/SceneOutlinerPanel.h"
#include "editor/panels/InspectorPanel.h"

// --- ECS 系统 ---
#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/function/render/RenderingSystem/RenderingSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h"
#include "runtime/function/physics/PhysicsScene.h"

// --- 资源管理 ---
#include "runtime/resource/resourceManager/resourceManager.h"
#include "runtime/function/res_type/Model/Model.h"
#include "runtime/function/res_type/shader/shader.h"
#include "runtime/function/res_type/Material/PBRMaterial.h"

namespace Lizeral {
    namespace Reflection {
        class TypeMetaRegister {
        public:
            static void Register(); 
        };
    }
}

// =====================================================================
// 【新增】：核心视口类，用于替代原来的黑屏 QLabel
// =====================================================================
class EngineViewportWidget : public QOpenGLWidget {
public:

    std::function<void()> onInitGL;
    
    EngineViewportWidget(Lizeral::Registry* registry, Lizeral::RenderingSystem* renderSys, QWidget* parent = nullptr) 
        : QOpenGLWidget(parent), m_Registry(registry), m_RenderSystem(renderSys) {


        // 允许此 Widget 接收键盘和鼠标事件
        setFocusPolicy(Qt::StrongFocus);
        setMouseTracking(true); 
    }

protected:
    // 1. OpenGL 上下文准备好时执行一次
    void initializeGL() override {
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

    // 2. 窗口大小改变时
    void resizeGL(int w, int h) override {
        glViewport(0, 0, w, h);
        // 未来可以在这里更新 CameraComponent 的 AspectRatio
    }

    // 3. 每帧渲染逻辑
    void paintGL() override {
        // 【关键修复】：获取 Qt 当前偷偷绑定的 FBO ID！
        // 你也可以用 this->defaultFramebufferObject()
        GLint qtCurrentFBO;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &qtCurrentFBO);
        
        // 告诉渲染系统：以后你的“屏幕”就是这个 ID！不要再绑 0 了！
        if (m_RenderSystem) {
            m_RenderSystem->SetDefaultFBO(static_cast<unsigned int>(qtCurrentFBO));
        }

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (m_RenderSystem && m_Registry) {
            m_RenderSystem->Tick(*m_Registry, 0.016f); 
        }
    }

    // ---------------------------------------------------------
    // 【预留接口】：将 Qt 的鼠标键盘事件转发给引擎 Input 系统
    // ---------------------------------------------------------
    void mousePressEvent(QMouseEvent *event) override {
        // 例：Lizeral::Input::GetInstance().SetMouseButtonDown(Left, true);
    }
    void mouseMoveEvent(QMouseEvent *event) override {
        // 例：Lizeral::Input::GetInstance().SetMousePosition(event->pos().x(), ...);
    }
    void keyPressEvent(QKeyEvent *event) override {
        // 例：Lizeral::Input::GetInstance().SetKeyDown(event->key(), true);
    }

private:
    Lizeral::Registry* m_Registry { nullptr };
    Lizeral::RenderingSystem* m_RenderSystem { nullptr };
};

// =====================================================================
// 编辑器主窗口
// =====================================================================
class LizeralEditorWindow : public QMainWindow {
public:
    LizeralEditorWindow() {
        setWindowTitle("Lizeral Engine - Full Architecture Test");
        resize(1366, 768);

        // 1. 初始化引擎反射系统与 UI 工厂
        Lizeral::Reflection::TypeMetaRegister::Register();
        Lizeral::DataTypeUIFactory::Initialize();

        // 2. 创建全局 ECS 注册表
        m_globalRegistry = new Lizeral::Registry();

        // 3. 构建编辑器 UI 架构
        setupUI();

        // // 4. 植入测试数据
        // populateTestData();

        // // 5. 初始化并启动引擎游戏循环
        // initEngineSystems();
    }

    ~LizeralEditorWindow() {
        delete m_globalRegistry;
    }

private:
    Lizeral::Registry* m_globalRegistry { nullptr };
    Lizeral::SceneOutlinerPanel* m_outlinerPanel { nullptr };
    Lizeral::InspectorPanel* m_inspectorPanel { nullptr };
    
    // 【新增】：持有视口指针
    EngineViewportWidget* m_viewportWidget { nullptr };

    // --- 引擎核心系统实例 ---
    Lizeral::PhysicsScene m_physicsScene;
    Lizeral::PhysicsSystem m_physicsSystem;
    Lizeral::RenderingSystem m_renderSystem;
    Lizeral::CameraSystem m_cameraSystem;
    
    QTimer* m_engineTimer { nullptr };
    QElapsedTimer m_timeTracker;
    float m_lastTime = 0.0f;

    void setupUI() {
        // --- A. 中央区域：真实的 OpenGL 视口 ---
        QWidget* centralWidget = new QWidget(this);
        centralWidget->setStyleSheet("background-color: #1e1e1e;"); 
        QVBoxLayout* centralLayout = new QVBoxLayout(centralWidget);
        centralLayout->setContentsMargins(0, 0, 0, 0); // 让 OpenGL 画面填充满中央区域

        // 创建真正的视口，并将 Registry 和 RenderSystem 的引用传给它
        m_viewportWidget = new EngineViewportWidget(m_globalRegistry, &m_renderSystem, centralWidget);

        m_viewportWidget->onInitGL = [this]() {
            std::cout << "[Editor] OpenGL Context Ready! Loading Assets..." << std::endl;
            this->populateTestData();   // 现在可以安全加载模型和 Shader 了
            this->initEngineSystems();  // 启动物理和渲染的 60FPS 循环
        };

        centralLayout->addWidget(m_viewportWidget);

        setCentralWidget(centralWidget);

        // --- B. 左侧区域：Scene Outliner ---
        QDockWidget* outlinerDock = new QDockWidget("Scene Outliner", this);
        outlinerDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
        
        m_outlinerPanel = new Lizeral::SceneOutlinerPanel(outlinerDock);
        m_outlinerPanel->SetRegistry(m_globalRegistry);
        
        outlinerDock->setWidget(m_outlinerPanel);
        addDockWidget(Qt::LeftDockWidgetArea, outlinerDock);

        // --- C. 右侧区域：Inspector ---
        QDockWidget* inspectorDock = new QDockWidget("Inspector", this);
        inspectorDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

        m_inspectorPanel = new Lizeral::InspectorPanel(inspectorDock);
        m_inspectorPanel->SetRegistry(m_globalRegistry);

        QScrollArea* scrollArea = new QScrollArea(inspectorDock);
        scrollArea->setWidgetResizable(true);
        scrollArea->setWidget(m_inspectorPanel);
        
        inspectorDock->setWidget(scrollArea);
        addDockWidget(Qt::RightDockWidgetArea, inspectorDock);

        // --- 信号绑定 ---
        connect(&Lizeral::EditorSelection::Get(), &Lizeral::EditorSelection::OnEntitySelected,
                m_inspectorPanel, &Lizeral::InspectorPanel::BindEntity);

        connect(&Lizeral::EditorSelection::Get(), &Lizeral::EditorSelection::OnEntityDataModified,
                this, [this](Lizeral::Entity) {
                    if (this->m_outlinerPanel) {
                        this->m_outlinerPanel->Refresh();
                    }
                });
    }

    void populateTestData() {
        // ==========================================
        // 1. 加载资源 (与你 Sandbox 中的逻辑一样)
        // ==========================================
        Lizeral::ResourceManager::GetInstance().SetRootPath("");
        
        auto boxModel = Lizeral::ResourceManager::GetInstance().Load<Lizeral::Model>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\box_with_uv.glb");
        auto pbrShader = std::make_shared<Lizeral::Shader>(
            "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.vert",
            "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.frag"
        );

        if (boxModel) {
            for (auto& mesh : boxModel->GetMeshes()) {
                auto pbrMat = std::dynamic_pointer_cast<Lizeral::PBRMaterial>(mesh.m_Material);
                if (pbrMat) {
                    pbrMat->SetShader(pbrShader);
                    pbrMat->m_Metallic = 0.1f;
                    pbrMat->m_Roughness = 0.9f;
                }
            }
        } else {
            std::cout << "[Editor Warning] Box model failed to load!" << std::endl;
        }

        // ==========================================
        // 2. 创建 ECS 实体
        // ==========================================
        // 测试实体 1：主摄像机
        // 1. 设置主摄像机
        Lizeral::Entity camEntity = m_globalRegistry->create();
        m_globalRegistry->emplace<Lizeral::NameComponent>(camEntity, "Main Camera");
        
        auto& camTrans = m_globalRegistry->emplace<Lizeral::TransformComponent>(camEntity);
        camTrans.setPosition(Lizeral::Vector3(0.0f, 2.0f, 15.0f)); 
        // 【关键修复】：由于没有 ControlSystem，我们直接硬编码让它看向 -Z 轴原点！
        camTrans.setForward(Lizeral::Vector3(0.0f, 0.0f, -1.0f)); 
        
        auto& cam = m_globalRegistry->emplace<Lizeral::CameraComponent>(camEntity);
        cam.setFov(45.0f);
        cam.setAspect(16.0f / 9.0f);
        cam.setzNear(0.1f);
        cam.setzFar(1000.0f);
        cam.setMain(true); 

        // 2. 设置光照（退回 Sandbox 里最稳妥的 setForward）
        Lizeral::Entity lightEntity = m_globalRegistry->create();
        m_globalRegistry->emplace<Lizeral::NameComponent>(lightEntity, "Directional Light");
        auto& lightTrans = m_globalRegistry->emplace<Lizeral::TransformComponent>(lightEntity);
        // 同样用 setForward 来决定光照方向
        lightTrans.setForward(Lizeral::Vector3(1.0f, 1.0f, 1.0f)); 
        m_globalRegistry->emplace<Lizeral::DirectionLightComponent>(lightEntity).setIntensity(3.0f);

        // 测试实体 3：放一个我们刚才加载的箱子！！！
        Lizeral::Entity boxEntity = m_globalRegistry->create();
        m_globalRegistry->emplace<Lizeral::NameComponent>(boxEntity, "Test Box");
        m_globalRegistry->emplace<Lizeral::TransformComponent>(boxEntity); // 默认在 (0,0,0)
        
        if (boxModel) {
            m_globalRegistry->emplace<Lizeral::ModelComponent>(boxEntity, boxModel);
        }

        // 刷新大纲以显示这些初始数据，并默认选中箱子
        m_outlinerPanel->Refresh();
        Lizeral::EditorSelection::Get().SelectEntity(boxEntity);
    }

    void initEngineSystems() {
        m_physicsScene.Initialize();
        m_physicsSystem.Initialize(&m_physicsScene);
        // 注意：m_renderSystem.Initialize() 已经在 EngineViewportWidget::initializeGL() 中调用了！

        // 设置游戏主循环 (约 60FPS = 16ms)
        m_engineTimer = new QTimer(this);
        // 在现代 Qt 中，普通成员函数可以直接用这种方式 connect，无需将其声明为 slots
        connect(m_engineTimer, &QTimer::timeout, this, &LizeralEditorWindow::EngineTick);
        
        m_timeTracker.start();
        m_engineTimer->start(16); 
    }

    void EngineTick() {
        float currentTime = m_timeTracker.elapsed() / 1000.0f;
        float deltaTime = currentTime - m_lastTime;
        m_lastTime = currentTime;
        if (deltaTime > 0.1f) deltaTime = 0.1f;

        // 执行引擎逻辑
        m_physicsSystem.Tick(deltaTime, *m_globalRegistry);
        m_cameraSystem.Tick(*m_globalRegistry);
        
        // 触发 OpenGL 视口的 paintGL
        if (m_viewportWidget) {
            m_viewportWidget->update();
        }
    }
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QPalette darkPalette;
    darkPalette.setColor(QPalette::Window, QColor(43, 43, 43));
    darkPalette.setColor(QPalette::WindowText, Qt::white);
    darkPalette.setColor(QPalette::Base, QColor(25, 25, 25));
    darkPalette.setColor(QPalette::AlternateBase, QColor(43, 43, 43));
    darkPalette.setColor(QPalette::ToolTipBase, Qt::white);
    darkPalette.setColor(QPalette::ToolTipText, Qt::white);
    darkPalette.setColor(QPalette::Text, Qt::white);
    darkPalette.setColor(QPalette::Button, QColor(43, 43, 43));
    darkPalette.setColor(QPalette::ButtonText, Qt::white);
    darkPalette.setColor(QPalette::BrightText, Qt::red);
    app.setPalette(darkPalette);

    LizeralEditorWindow editorWindow;
    editorWindow.show();

    return app.exec();
}