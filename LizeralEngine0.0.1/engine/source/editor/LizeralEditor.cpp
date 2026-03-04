#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <QApplication>
#include <QMainWindow>
#include <QDockWidget>
#include <QScrollArea>
#include <QOpenGLWidget>
#include <QVBoxLayout>
#include <QLabel>
#include <iostream>

// --- 引擎侧模块 ---
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/componentAll.h"

// --- 编辑器侧模块 ---
#include "editor/factory/DataTypeUIFactory/DataTypeFactory.h"
#include "editor/selection/EditorSelection.h"
#include "editor/panels/SceneOutlinerPanel.h"
#include "editor/panels/InspectorPanel.h"

namespace Lizeral {
    namespace Reflection {
        class TypeMetaRegister {
        public:
            static void Register(); 
        };
    }
}

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

        // 4. 植入测试数据
        populateTestData();
    }

    ~LizeralEditorWindow() {
        delete m_globalRegistry;
    }

private:
    Lizeral::Registry* m_globalRegistry { nullptr };
    Lizeral::SceneOutlinerPanel* m_outlinerPanel { nullptr };
    Lizeral::InspectorPanel* m_inspectorPanel { nullptr };

    void setupUI() {
        // --- A. 中央区域：占位视口 ---
        QWidget* centralWidget = new QWidget(this);
        centralWidget->setStyleSheet("background-color: #1e1e1e;"); 
        QVBoxLayout* centralLayout = new QVBoxLayout(centralWidget);
        QLabel* viewportLabel = new QLabel("Engine Viewport\n(OpenGL Rendering Area)", centralWidget);
        viewportLabel->setAlignment(Qt::AlignCenter);
        viewportLabel->setStyleSheet("color: #666666; font-size: 24px; font-weight: bold;");
        centralLayout->addWidget(viewportLabel);
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

        // ==========================================================
        // 核心信号绑定：打通任督二脉
        // SceneOutliner -> (点击) -> EditorSelection -> (发射) -> Inspector
        // ==========================================================
        connect(&Lizeral::EditorSelection::Get(), &Lizeral::EditorSelection::OnEntitySelected,
                m_inspectorPanel, &Lizeral::InspectorPanel::BindEntity);
    }

    void populateTestData() {
        // 测试实体 1：隐身摄像机（测试 EditorOnlyComponent 的黑名单过滤逻辑）
        Lizeral::Entity camEntity = m_globalRegistry->create();
        m_globalRegistry->emplace<Lizeral::NameComponent>(camEntity, "Editor Camera");
        m_globalRegistry->emplace<Lizeral::EditorOnlyComponent>(camEntity); // 应该被大纲过滤掉！
        m_globalRegistry->emplace<Lizeral::CameraComponent>(camEntity);
        m_globalRegistry->emplace<Lizeral::TransformComponent>(camEntity);

        // 测试实体 2：普通光照
        Lizeral::Entity lightEntity = m_globalRegistry->create();
        m_globalRegistry->emplace<Lizeral::NameComponent>(lightEntity, "Directional Light");
        m_globalRegistry->emplace<Lizeral::TransformComponent>(lightEntity);
        m_globalRegistry->emplace<Lizeral::DirectionLightComponent>(lightEntity);

        // 测试实体 3：物理盒子
        Lizeral::Entity cubeEntity = m_globalRegistry->create();
        m_globalRegistry->emplace<Lizeral::NameComponent>(cubeEntity, "Physics Cube");
        m_globalRegistry->emplace<Lizeral::TransformComponent>(cubeEntity);
        m_globalRegistry->emplace<Lizeral::RigidBodyComponent>(cubeEntity);
        m_globalRegistry->emplace<Lizeral::ColliderComponent>(cubeEntity);

        // 刷新大纲以显示这些初始数据
        m_outlinerPanel->Refresh();

        // 默认选中第一个可视实体（光照）
        Lizeral::EditorSelection::Get().SelectEntity(lightEntity);
    }
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    // 设置全局极暗主题（可选，让编辑器看起来更专业）
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