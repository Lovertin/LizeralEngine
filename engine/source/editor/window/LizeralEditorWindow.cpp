#include "LizeralEditorWindow.h"
#include <QIcon>
#include <QShortcut>
#include <QFileDialog>
#include <QMessageBox>
#include <QFile>
#include <QTextStream>
#include <iostream>
#include <QFileSystemWatcher>

LizeralEditorWindow::LizeralEditorWindow(){

    setWindowTitle("Lizeral Engine");
    resize(1366, 768);

    setWindowIcon(QIcon("C:/Lizeral Engine/LizeralEngine0.0.1/asset/logos/Lizeral_Engine.jpg"));

    Lizeral::Reflection::TypeMetaRegister::Register();
    Lizeral::DataTypeUIFactory::Initialize();

    m_globalRegistry = new Lizeral::Registry();

    Lizeral::EditorContext::Get().SetRegistry(m_globalRegistry);

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

void LizeralEditorWindow::loadStyleSheet(const QString& sheetName) {
    QFile file(sheetName);
    if (file.open(QFile::ReadOnly | QFile::Text)) {
        QTextStream stream(&file);
        QString qss = stream.readAll();
        this->setStyleSheet(qss);
        file.close();
        std::cout << "[Editor] Loaded Style Sheet: " << sheetName.toStdString() << std::endl;
    } else {
        std::cerr << "[Editor] Failed to load Style Sheet: " << sheetName.toStdString() << std::endl;
    }
}

void LizeralEditorWindow::setupUI(){
    m_controlPanel = new Lizeral::EditorControlPanel(this);
    addToolBar(Qt::TopToolBarArea, m_controlPanel);

    loadStyleSheet("C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/editor/style/UE5_Dark.qss");

    QString qssPath = "C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/editor/style/UE5_Dark.qss";
    loadStyleSheet(qssPath);

    QFileSystemWatcher* watcher = new QFileSystemWatcher(this);
    watcher->addPath(qssPath);
    
    connect(watcher, &QFileSystemWatcher::fileChanged, this, [this, watcher, qssPath](const QString& path) {
        std::cout << "[Editor] Hot-Reloading Style Sheet..." << std::endl;
        loadStyleSheet(path);
        if (!watcher->files().contains(qssPath)) {
            watcher->addPath(qssPath);
        }
    });

    // Listen for the Save signal
    connect(m_controlPanel, &Lizeral::EditorControlPanel::OnSaveScene, this, [this]() {
        QString filePath = QFileDialog::getSaveFileName(
            this, "Save Scene", "", "JSON Files (*.json)");
            
        if (!filePath.isEmpty()) {
            Lizeral::SceneSerializer::SaveToFile(m_globalRegistry, filePath.toStdString());
        }
    });

    connect(m_controlPanel, &Lizeral::EditorControlPanel::OnLoadScene, this, [this]() {
        QString filePath = QFileDialog::getOpenFileName(
            this, "Load Scene", "", "JSON Files (*.json)");
            
        if (!filePath.isEmpty()) {
            m_renderSystem.WaitIdle();
            // load Scene
            Lizeral::SceneSerializer::LoadFromFile(filePath.toStdString(), m_globalRegistry);
            
            Lizeral::EditorContext::Get().GetCommandManager()->Clear();
            if (m_outlinerPanel) m_outlinerPanel->Refresh();
            Lizeral::EditorSelection::Get().SelectEntity(Lizeral::null_entity);
        }
    });

    connect(m_controlPanel, &Lizeral::EditorControlPanel::OnPlayModeEntered, this, [this]() {
        std::cout << "[Editor] Saving Scene Snapshot..." << std::endl;

        m_playModeSnapshot = Lizeral::SceneSerializer::Serialize(m_globalRegistry);
    });

    connect(m_controlPanel, &Lizeral::EditorControlPanel::OnEditModeEntered, this, [this]() {
        std::cout << "[Editor] Restoring Scene Snapshot..." << std::endl;

        m_renderSystem.WaitIdle();

        m_physicsSystem.ResetPhysicsWorld();
        Lizeral::SceneSerializer::Deserialize(m_playModeSnapshot, m_globalRegistry);

        if (m_outlinerPanel) m_outlinerPanel->Refresh();
        Lizeral::EditorSelection::Get().SelectEntity(Lizeral::null_entity); 
    });

    QWidget* centralWidget = new QWidget(this);
    centralWidget->setStyleSheet("background-color: #1e1e1e;"); 
    QVBoxLayout* centralLayout = new QVBoxLayout(centralWidget);
    centralLayout->setContentsMargins(0, 0, 0, 0);
    m_renderSystem.SetRenderPipelinePreset(Lizeral::VulkanRenderingSystem::RenderPipelinePreset::Stable);

    m_viewportWidget = new Lizeral::EngineViewportWidget(m_globalRegistry, &m_renderSystem, centralWidget);

    m_viewportWidget->SetPhysicsSystem(&m_physicsSystem);

    m_viewportWidget->onInitVulkan = [this]() {
        std::cout << "[Editor] Vulkan Context Ready! Loading Assets..." << std::endl;
        this->populateTestData();  
        this->initEngineSystems();  
    };

    centralLayout->addWidget(m_viewportWidget);

    setCentralWidget(centralWidget);

    // --- Scene Outliner ---
    QDockWidget* outlinerDock = new QDockWidget("Scene Outliner", this);
    outlinerDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    
    m_outlinerPanel = new Lizeral::SceneOutlinerPanel(outlinerDock);
    m_outlinerPanel->SetRegistry(m_globalRegistry);
    
    outlinerDock->setWidget(m_outlinerPanel);
    addDockWidget(Qt::LeftDockWidgetArea, outlinerDock);

    // --- Inspector ---
    QDockWidget* inspectorDock = new QDockWidget("Inspector", this);
    inspectorDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    m_inspectorPanel = new Lizeral::InspectorPanel(inspectorDock);
    m_inspectorPanel->SetRegistry(m_globalRegistry);
    m_inspectorPanel->setMinimumWidth(320);

    QScrollArea* scrollArea = new QScrollArea(inspectorDock);
    scrollArea->setWidgetResizable(true);
    scrollArea->setWidget(m_inspectorPanel);
    
    inspectorDock->setWidget(scrollArea);
    addDockWidget(Qt::RightDockWidgetArea, inspectorDock);

    connect(&Lizeral::EditorSelection::Get(), &Lizeral::EditorSelection::OnEntitySelected,
            m_inspectorPanel, &Lizeral::InspectorPanel::BindEntity);

    connect(&Lizeral::EditorSelection::Get(), &Lizeral::EditorSelection::OnSubMeshSelectionChanged,
            this, [this](Lizeral::Entity entity, int32_t) {
                if (this->m_inspectorPanel) {
                    this->m_inspectorPanel->BindEntity(entity);
                }
            });

    connect(&Lizeral::EditorSelection::Get(), &Lizeral::EditorSelection::OnEntityDataModified,
            this, [this](Lizeral::Entity) {
                if (this->m_outlinerPanel) {
                    this->m_outlinerPanel->Refresh();
                }
            });

    QShortcut* undoShortcut = new QShortcut(QKeySequence("Ctrl+Z"), this);

    undoShortcut->setContext(Qt::ApplicationShortcut); 
    
    connect(undoShortcut, &QShortcut::activated, this, [this]() {
        auto cmdMgr = Lizeral::EditorContext::Get().GetCommandManager();
        if (cmdMgr->GetUndoSize() > 0) {
            std::cout << "[Command] Undo Execution... Remaining Stack: " << cmdMgr->GetUndoSize() - 1 << std::endl;
            cmdMgr->Undo();

            Lizeral::Entity current = Lizeral::EditorSelection::Get().GetSelected();
            if (current != Lizeral::null_entity) {
                Lizeral::EditorSelection::Get().SelectEntity(Lizeral::null_entity); 
                Lizeral::EditorSelection::Get().SelectEntity(current); 
            }
        } else {
            std::cout << "[Command] Undo Stack is Empty!" << std::endl;
        }
    });

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

    resizeDocks({outlinerDock, inspectorDock}, {250, 320}, Qt::Horizontal);
}

void LizeralEditorWindow::populateTestData(){

    // Lizeral::ResourceManager::GetInstance().SetRootPath("");

    Lizeral::Entity camEntity = m_globalRegistry->create();
    m_globalRegistry->emplace<Lizeral::NameComponent>(camEntity, "Camera");
    auto& camTrans = m_globalRegistry->emplace<Lizeral::TransformComponent>(camEntity);
    camTrans.setPosition(Lizeral::Vector3(0.0f, 2.0f, 15.0f)); 
    camTrans.setForward(Lizeral::Vector3(0.0f, 0.0f, -1.0f)); 
    auto& cam = m_globalRegistry->emplace<Lizeral::CameraComponent>(camEntity);
    cam.setFov(45.0f); cam.setAspect(16.0f / 9.0f); cam.setzNear(0.1f); cam.setzFar(1000.0f); cam.setMain(true); 
    auto& camCtrl = m_globalRegistry->emplace<Lizeral::CameraControlComponent>(camEntity);
    camCtrl.setYaw(0.0f); camCtrl.setPitch(0.0f); 
    camCtrl.setMoveSpeed(5.0f);

    m_outlinerPanel->Refresh();
    Lizeral::EditorSelection::Get().SelectEntity(camEntity);
}

void LizeralEditorWindow::initEngineSystems(){

    m_physicsScene.Initialize();
        m_physicsSystem.Initialize(&m_physicsScene);
        m_engineTimer = new QTimer(this);
        connect(m_engineTimer, &QTimer::timeout, this, &LizeralEditorWindow::EngineTick);
        
        m_timeTracker.start();
        m_engineTimer->start(8); //define fps
        //start(16) ~ 60fps     start(8) ~ 120fps
}

void LizeralEditorWindow::EngineTick()
{
    float currentTime = m_timeTracker.elapsed() / 1000.0f;
    float deltaTime = currentTime - m_lastTime;
    m_lastTime = currentTime;
    if (deltaTime > 0.1f) deltaTime = 0.1f;

    if(Lizeral::EditorContext::Get().getMode() == Lizeral::EditorMode::Play){
        m_physicsSystem.Tick(deltaTime, *m_globalRegistry);
    }

    m_physicsSystem.UpdateDebugLines(*m_globalRegistry);

    m_cameraControlSystem.Tick(deltaTime,*m_globalRegistry);
    m_cameraSystem.Tick(*m_globalRegistry);

    if (m_viewportWidget) {
        m_viewportWidget->SetDebugLines(m_physicsSystem.GetDebugLines());
        m_viewportWidget->update();
    }

    Lizeral::Input::GetInstance().Tick();

}
