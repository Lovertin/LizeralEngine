#include "LizeralEditorWindow.h"
#include <QShortcut>

LizeralEditorWindow::LizeralEditorWindow(){

    setWindowTitle("Lizeral Engine - Full Architecture Test");
    resize(1366, 768);

    // 1. 初始化引擎反射系统与 UI 工厂
    Lizeral::Reflection::TypeMetaRegister::Register();
    Lizeral::DataTypeUIFactory::Initialize();

    // 2. 创建全局 ECS 注册表
    m_globalRegistry = new Lizeral::Registry();

    Lizeral::EditorContext::Get().SetRegistry(m_globalRegistry);

    // 3. 构建编辑器 UI 架构
    setupUI();

}

LizeralEditorWindow::~LizeralEditorWindow(){
    if (m_engineTimer) {
        m_engineTimer->stop();
    }

    if (m_viewportWidget) {
        delete m_viewportWidget; 
        m_viewportWidget = nullptr;
    }

    if (m_globalRegistry) {
        delete m_globalRegistry;
        m_globalRegistry = nullptr;
    }
}

void LizeralEditorWindow::setupUI(){

    // --- A. 中央区域：真实的 OpenGL 视口 ---
    QWidget* centralWidget = new QWidget(this);
    centralWidget->setStyleSheet("background-color: #1e1e1e;"); 
    QVBoxLayout* centralLayout = new QVBoxLayout(centralWidget);
    centralLayout->setContentsMargins(0, 0, 0, 0); // 让 OpenGL 画面填充满中央区域

    // 创建真正的视口，并将 Registry 和 RenderSystem 的引用传给它
    m_viewportWidget = new Lizeral::EngineViewportWidget(m_globalRegistry, &m_renderSystem, centralWidget);

    m_viewportWidget->SetPhysicsSystem(&m_physicsSystem);

    m_viewportWidget->onInitVulkan = [this]() {
        std::cout << "[Editor] Vulkan Context Ready! Loading Assets..." << std::endl;
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

    // 1. 全局撤销快捷键
    QShortcut* undoShortcut = new QShortcut(QKeySequence("Ctrl+Z"), this);
    // 【核心修复1】：设置为 ApplicationShortcut，强行抢夺 QLineEdit 的 Ctrl+Z 拦截权！
    undoShortcut->setContext(Qt::ApplicationShortcut); 
    
    connect(undoShortcut, &QShortcut::activated, this, [this]() {
        auto cmdMgr = Lizeral::EditorContext::Get().GetCommandManager();
        if (cmdMgr->GetUndoSize() > 0) {
            std::cout << "[Command] Undo Execution... Remaining Stack: " << cmdMgr->GetUndoSize() - 1 << std::endl;
            cmdMgr->Undo();

            // 【核心修复2】：强制刷新 Inspector 面板（绕过 m_CurrentEntity 的拦截）
            Lizeral::Entity current = Lizeral::EditorSelection::Get().GetSelected();
            if (current != Lizeral::null_entity) {
                Lizeral::EditorSelection::Get().SelectEntity(Lizeral::null_entity); // 先置空
                Lizeral::EditorSelection::Get().SelectEntity(current);              // 再重绑
            }
        } else {
            std::cout << "[Command] Undo Stack is Empty!" << std::endl;
        }
    });

    // 2. 全局重做快捷键
    QShortcut* redoShortcut = new QShortcut(QKeySequence("Ctrl+Y"), this);
    redoShortcut->setContext(Qt::ApplicationShortcut);
    
    connect(redoShortcut, &QShortcut::activated, this, [this]() {
        auto cmdMgr = Lizeral::EditorContext::Get().GetCommandManager();
        if (cmdMgr->GetRedoSize() > 0) {
            std::cout << "[Command] Redo Execution..." << std::endl;
            cmdMgr->Redo();

            Lizeral::Entity current = Lizeral::EditorSelection::Get().GetSelected();
            if (current != Lizeral::null_entity) {
                Lizeral::EditorSelection::Get().SelectEntity(Lizeral::null_entity);
                Lizeral::EditorSelection::Get().SelectEntity(current);
            }
        }
    });
}

void LizeralEditorWindow::populateTestData(){

    Lizeral::ResourceManager::GetInstance().SetRootPath("");

    Lizeral::Entity camEntity = m_globalRegistry->create();
    m_globalRegistry->emplace<Lizeral::NameComponent>(camEntity, "Main Camera");
    auto& camTrans = m_globalRegistry->emplace<Lizeral::TransformComponent>(camEntity);
    camTrans.setPosition(Lizeral::Vector3(0.0f, 2.0f, 15.0f)); 
    camTrans.setForward(Lizeral::Vector3(0.0f, 0.0f, -1.0f)); 
    auto& cam = m_globalRegistry->emplace<Lizeral::CameraComponent>(camEntity);
    cam.setFov(45.0f); cam.setAspect(16.0f / 9.0f); cam.setzNear(0.1f); cam.setzFar(1000.0f); cam.setMain(true); 
    auto& camCtrl = m_globalRegistry->emplace<Lizeral::CameraControlComponent>(camEntity);
    camCtrl.setYaw(0.0f); camCtrl.setPitch(0.0f); 
    camCtrl.setMoveSpeed(5.0f);

    Lizeral::Entity lightEntity = m_globalRegistry->create();
    m_globalRegistry->emplace<Lizeral::NameComponent>(lightEntity, "Directional Light");
    auto& lightTrans = m_globalRegistry->emplace<Lizeral::TransformComponent>(lightEntity);
    lightTrans.setForward(Lizeral::Vector3(1.0f, 1.0f, 1.0f)); 
    auto& lighttype=m_globalRegistry->emplace<Lizeral::DirectionLightComponent>(lightEntity);
    lighttype.setIntensity(3.0f);

    Lizeral::Entity boxEntity = m_globalRegistry->create();
    m_globalRegistry->emplace<Lizeral::NameComponent>(boxEntity, "Test Box");
    m_globalRegistry->emplace<Lizeral::TransformComponent>(boxEntity); 

    Lizeral::Entity roomEntity = m_globalRegistry->create();
    m_globalRegistry->emplace<Lizeral::NameComponent>(roomEntity,"classroom");
    auto& roomTrans = m_globalRegistry->emplace<Lizeral::TransformComponent>(roomEntity);
    roomTrans.setPosition(Lizeral::Vector3(0.0f, -1.0f, 10.0f));
    roomTrans.setScale(Lizeral::Vector3(0.05f, 0.05f, 0.05f)); 
    m_globalRegistry->emplace<Lizeral::VulkanModelComponent>(roomEntity).setVulkanModelPath("C:/Lizeral Engine/LizeralEngine0.0.1/asset/scene_without_window.glb"); 

    Lizeral::Entity maserati = m_globalRegistry->create();
    m_globalRegistry->emplace<Lizeral::NameComponent>(maserati,"Maserati");
    auto& carTrans = m_globalRegistry->emplace<Lizeral::TransformComponent>(maserati);
    carTrans.setPosition(Lizeral::Vector3(80.0f,0.0f,0.0f));
    carTrans.setScale(Lizeral::Vector3(1.0f,1.0f,1.0f));
    m_globalRegistry->emplace<Lizeral::VulkanModelComponent>(maserati).setVulkanModelPath("C:/Lizeral Engine/LizeralEngine0.0.1/asset/maserati_withoutwindow.glb");

    // 刷新大纲以显示这些初始数据，并默认选中箱子
    m_outlinerPanel->Refresh();
    Lizeral::EditorSelection::Get().SelectEntity(boxEntity);
}

void LizeralEditorWindow::initEngineSystems(){

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

void LizeralEditorWindow::EngineTick()
{
    float currentTime = m_timeTracker.elapsed() / 1000.0f;
    float deltaTime = currentTime - m_lastTime;
    m_lastTime = currentTime;
    if (deltaTime > 0.1f) deltaTime = 0.1f;

    // 执行引擎逻辑
    m_physicsSystem.Tick(deltaTime, *m_globalRegistry);
    m_cameraControlSystem.Tick(deltaTime,*m_globalRegistry);
    m_cameraSystem.Tick(*m_globalRegistry);

    if (m_viewportWidget) {
        m_viewportWidget->update();
    }

    Lizeral::Input::GetInstance().Tick();

}