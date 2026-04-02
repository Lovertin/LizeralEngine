#pragma once

#include <glad/glad.h>
// 注意：在 Qt 中，我们尽量不使用 GLFW，因为窗口和输入由 Qt 接管了
// 但如果你的底层代码(比如 Input)还强依赖它，可以保留。
#include <GLFW/glfw3.h>

#include <QApplication>
#include <QMainWindow>
#include <QDockWidget>
#include <QScrollArea>
#include <QOpenGLWidget>
#include <QOpenGLContext>
#include <QVBoxLayout>
#include <QLabel>
#include <iostream>
#include <QTimer>
#include <QElapsedTimer>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QCursor>

#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/componentAll.h"

#include "factory/DataTypeUIFactory/DataTypeFactory.h"
#include "selection/EditorSelection.h"
#include "panels/SceneOutlinerPanel.h"
#include "panels/InspectorPanel.h"

#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h"
#include "runtime/function/physics/PhysicsScene.h"
#include "runtime/function/render/CameraControlSystem/CameraControlSystem.h"
#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystem.h"
#include "runtime/function/input/input.h"

#include "runtime/resource/resourceManager/resourceManager.h"