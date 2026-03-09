#include "LizeralEditorWindow.h"

LizeralEditorWindow::LizeralEditorWindow(){

    setWindowTitle("Lizeral Engine - Full Architecture Test");
    resize(1366, 768);

    // 1. 初始化引擎反射系统与 UI 工厂
    Lizeral::Reflection::TypeMetaRegister::Register();
    Lizeral::DataTypeUIFactory::Initialize();

    // 2. 创建全局 ECS 注册表
    m_globalRegistry = new Lizeral::Registry();

    // 3. 构建编辑器 UI 架构
    setupUI();

}

void LizeralEditorWindow::setupUI(){

    // --- A. 中央区域：真实的 OpenGL 视口 ---
    QWidget* centralWidget = new QWidget(this);
    centralWidget->setStyleSheet("background-color: #1e1e1e;"); 
    QVBoxLayout* centralLayout = new QVBoxLayout(centralWidget);
    centralLayout->setContentsMargins(0, 0, 0, 0); // 让 OpenGL 画面填充满中央区域

    // 创建真正的视口，并将 Registry 和 RenderSystem 的引用传给它
    m_viewportWidget = new EngineViewportWidget(m_globalRegistry, &m_renderSystem, centralWidget);

    m_viewportWidget->SetPhysicsSystem(&m_physicsSystem);

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

void LizeralEditorWindow::populateTestData(){

    Lizeral::ResourceManager::GetInstance().SetRootPath("");
        
    // 提前准备好一个公用的 Shader
    auto pbrShader = std::make_shared<Lizeral::Shader>(
        "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.vert",
        "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.frag"
    );

    // ==========================================
    // 1. 创建摄像机和灯光 (保持你原来的代码不变)
    // ==========================================
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
    
    // A. 挂载 ModelComponent (调用默认构造函数)
    auto& modelComp = m_globalRegistry->emplace<Lizeral::ModelComponent>(boxEntity);
    
    // B. 设置路径并让组件自己去加载
    modelComp.m_ModelPath = "C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\monkey.glb";
    modelComp.LoadResources(); // 这会自动加载模型，并准备好 override 数组！

    // C. 体验材质覆盖！定制这个独一无二的箱子
    if (modelComp.m_Model) {
        for (auto& mesh : modelComp.m_Model->GetMeshes()) {
            auto pbrMat = std::dynamic_pointer_cast<Lizeral::PBRMaterial>(mesh.m_Material);
            if (pbrMat) {
                // 仅仅赋予 Shader，绝不修改它从 Blender 里继承来的 m_Albedo 颜色！
                pbrMat->SetShader(pbrShader);
            }
        }
    }

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
    
    // 触发 OpenGL 视口的 paintGL
    if (m_viewportWidget) {
        m_viewportWidget->update();
    }

    Lizeral::Input::GetInstance().Tick();

}