#pragma once 
#include "editor/EditorHeader.h"
#include "editor/viewport/EngineViewportWidget.h"
#include "editor/context/EditorContext.h"
#include "editor/panels/EditorControlPanel.h"

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
    LizeralEditorWindow();

    ~LizeralEditorWindow();

private:
    Lizeral::Registry* m_globalRegistry { nullptr };
    Lizeral::SceneOutlinerPanel* m_outlinerPanel { nullptr };
    Lizeral::InspectorPanel* m_inspectorPanel { nullptr };
    Lizeral::EditorControlPanel* m_controlPanel { nullptr };
    
    // 【新增】：持有视口指针
    Lizeral::EngineViewportWidget* m_viewportWidget { nullptr };

    // --- 引擎核心系统实例 ---
    Lizeral::PhysicsScene m_physicsScene;
    Lizeral::PhysicsSystem m_physicsSystem;
    Lizeral::VulkanRenderingSystem m_renderSystem;
    Lizeral::CameraSystem m_cameraSystem;
    Lizeral::CameraControlSystem m_cameraControlSystem; 
    
    QTimer* m_engineTimer { nullptr };
    QElapsedTimer m_timeTracker;
    float m_lastTime = 0.0f;

    void setupUI();

    void populateTestData();

    void initEngineSystems();

    void EngineTick();
};