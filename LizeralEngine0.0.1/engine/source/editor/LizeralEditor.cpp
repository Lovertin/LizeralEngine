// #include <glad/glad.h>
// // 注意：在 Qt 中，我们尽量不使用 GLFW，因为窗口和输入由 Qt 接管了
// // 但如果你的底层代码(比如 Input)还强依赖它，可以保留。
// #include <GLFW/glfw3.h>

// #include <QApplication>
// #include <QMainWindow>
// #include <QDockWidget>
// #include <QScrollArea>
// #include <QOpenGLWidget>
// #include <QOpenGLContext> // 【新增】：用于获取 Qt 的 OpenGL 函数指针
// #include <QVBoxLayout>
// #include <QLabel>
// #include <iostream>
// #include <QTimer>
// #include <QElapsedTimer>
// #include <QMouseEvent>
// #include <QKeyEvent>

// // --- 引擎侧模块 ---
// #include "runtime/core/meta/reflection/reflection.h"
// #include "runtime/function/ecs/registry.h"
// #include "runtime/function/ecs/components/componentAll.h"

// // --- 编辑器侧模块 ---
// #include "editor/factory/DataTypeUIFactory/DataTypeFactory.h"
// #include "editor/selection/EditorSelection.h"
// #include "editor/panels/SceneOutlinerPanel.h"
// #include "editor/panels/InspectorPanel.h"

// // --- ECS 系统 ---
// #include "runtime/function/physics/PhysicsSystem.h"
// #include "runtime/function/render/RenderingSystem/RenderingSystem.h"
// #include "runtime/function/render/CameraSystem/CameraSystem.h"
// #include "runtime/function/physics/PhysicsScene.h"

// // --- 资源管理 ---
// #include "runtime/resource/resourceManager/resourceManager.h"
// #include "runtime/function/res_type/Model/Model.h"
// #include "runtime/function/res_type/shader/shader.h"
// #include "runtime/function/res_type/Material/PBRMaterial.h"

#include "EditorHeader.h"
#include "viewport/EngineViewportWidget.h"
#include "window/LizeralEditorWindow.h"


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
