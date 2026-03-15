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

#include "LizeralEditorWindow.h"
#include <QResizeEvent>

// ... 构造函数保持不变 ...

void LizeralEditorWindow::setupUI(){

    // ==========================================================
    // Layer 0 (底层)：真正的 OpenGL/Vulkan 视口，永远铺满整个主窗口
    // ==========================================================
    m_viewportWidget = new Lizeral::EngineViewportWidget(m_globalRegistry, &m_renderSystem, this);
    
    m_viewportWidget->onInitVulkan = [this]() {
        std::cout << "[Editor] Vulkan Context Ready! Loading Assets..." << std::endl;
        this->populateTestData();  
        this->initEngineSystems();  
    };

    // 直接将视口设为唯一的中央组件
    setCentralWidget(m_viewportWidget);


    // ==========================================================
    // Layer 1 (顶层)：悬浮在视口上的 UI 面板 (作为 viewport 的子节点)
    // ==========================================================

    // --- 左侧区域：Scene Outliner ---
    m_outlinerPanel = new Lizeral::SceneOutlinerPanel(m_viewportWidget); // 注意：父节点传 viewport
    m_outlinerPanel->SetRegistry(m_globalRegistry);
    // 必须设置不透明的背景色，否则文字会和 3D 场景糊在一起
    m_outlinerPanel->setStyleSheet("background-color: rgba(30, 30, 30, 240); border: 1px solid #444;");
    m_outlinerPanel->show(); // 强制显示子窗口

    // --- 右侧区域：Inspector ---
    m_inspectorPanel = new Lizeral::InspectorPanel(m_viewportWidget);
    m_inspectorPanel->SetRegistry(m_globalRegistry);
    m_inspectorPanel->setStyleSheet("background-color: transparent;");

    m_inspectorScrollArea = new QScrollArea(m_viewportWidget); // 注意：父节点传 viewport
    m_inspectorScrollArea->setWidgetResizable(true);
    m_inspectorScrollArea->setWidget(m_inspectorPanel);
    m_inspectorScrollArea->setStyleSheet("background-color: rgba(30, 30, 30, 240); border: 1px solid #444;");
    m_inspectorScrollArea->show(); // 强制显示

    // --- 信号绑定 (完全不用动) ---
    connect(&Lizeral::EditorSelection::Get(), &Lizeral::EditorSelection::OnEntitySelected,
            m_inspectorPanel, &Lizeral::InspectorPanel::BindEntity);

    connect(&Lizeral::EditorSelection::Get(), &Lizeral::EditorSelection::OnEntityDataModified,
            this, [this](Lizeral::Entity) {
                if (this->m_outlinerPanel) {
                    this->m_outlinerPanel->Refresh();
                }
            });
}

// ==========================================================
// 核心排版逻辑：手动将悬浮面板钉在左右两侧
// ==========================================================
void LizeralEditorWindow::resizeEvent(QResizeEvent* event) {
    QMainWindow::resizeEvent(event);
    
    // 如果视口还没建好，直接跳过
    if (!m_viewportWidget) return;

    int panelWidth = 320; // 你想要的面板宽度
    int padding = 0;      // 面板距离窗口边缘的空隙（如果想完全贴边可以设为0）

    // 钉住左边 Outliner
    if (m_outlinerPanel) {
        m_outlinerPanel->setGeometry(padding, padding, 
                                     panelWidth, this->height() - padding * 2);
    }

    // 钉住右边 Inspector
    if (m_inspectorScrollArea) {
        m_inspectorScrollArea->setGeometry(this->width() - panelWidth - padding, padding, 
                                           panelWidth, this->height() - padding * 2);
    }
}

void LizeralEditorWindow::populateTestData(){

    Lizeral::ResourceManager::GetInstance().SetRootPath("");

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
    m_globalRegistry->emplace<Lizeral::NameComponent>(boxEntity, "Test Box (Vulkan)");
    
    auto& boxTrans = m_globalRegistry->emplace<Lizeral::TransformComponent>(boxEntity); 
    boxTrans.setPosition(Lizeral::Vector3(0.0f, 0.0f, 0.0f));
    
    // ★ 挂载 Vulkan 专用的模型组件
    auto& vkModelComp = m_globalRegistry->emplace<Lizeral::VulkanModelComponent>(boxEntity);
    
    // 传入你的测试模型路径 (MeshletBuilder 会自动在底层接管并解析它的材质)
    vkModelComp.setVulkanModelPath("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\scene_without_window.glb"); 
    // 或者用你的猴头："C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\monkey.glb"

    // 刷新大纲以显示这些初始数据，并默认选中它
    m_outlinerPanel->Refresh();
    Lizeral::EditorSelection::Get().SelectEntity(boxEntity);
}

void LizeralEditorWindow::initEngineSystems(){

    m_physicsScene.Initialize();
    m_physicsSystem.Initialize(&m_physicsScene);
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