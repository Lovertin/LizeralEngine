#include <iostream>
#include <exception>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

// --- Vulkan & Engine RHI ---
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"

// --- ★ 引入我们新封装的渲染系统 ★ ---
#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystem.h"

// --- ECS & Camera System ---
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/ecs/components/Model/VulkanModelComponent.h" 
#include "runtime/function/input/input.h"
#include "runtime/function/render/CameraControlSystem/CameraControlSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h" 
#include "runtime/core/math/matrix4.h"

using namespace Lizeral;

const uint32_t WIDTH = 1280;
const uint32_t HEIGHT = 720;

// ====================================================================
// 工具函数与回调区 (保持不变)
// ====================================================================
void glfwKeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    Key lizKey; bool found = true;
    switch (key) {
        case GLFW_KEY_W: lizKey = Key::W; break;
        case GLFW_KEY_A: lizKey = Key::A; break;
        case GLFW_KEY_S: lizKey = Key::S; break;
        case GLFW_KEY_D: lizKey = Key::D; break;
        case GLFW_KEY_Q: lizKey = Key::Q; break;
        case GLFW_KEY_E: lizKey = Key::E; break;
        case GLFW_KEY_LEFT_SHIFT: lizKey = Key::LEFT_SHIFT; break;
        case GLFW_KEY_ESCAPE: lizKey = Key::ESC; break;
        default: found = false; break;
    }
    if (found) Input::GetInstance().SetKeyDown(lizKey, action != GLFW_RELEASE);
}

void glfwCursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    Input::GetInstance().SetMousePosition(static_cast<float>(xpos), static_cast<float>(ypos));
}

void glfwMouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        Input::GetInstance().SetMouseButtonDown(MouseButton::Right, action != GLFW_RELEASE);
    }
}

// ====================================================================
// 主程序入口
// ====================================================================
int main() {
    GLFWwindow* window = nullptr;
    VkSurfaceKHR surface = VK_NULL_HANDLE;
    VkInstance instance = VK_NULL_HANDLE;

    try {
        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API); 
        glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
        window = glfwCreateWindow(WIDTH, HEIGHT, "Lizeral Engine - Refactored Vulkan Rendering System", nullptr, nullptr);

        glfwSetKeyCallback(window, glfwKeyCallback);
        glfwSetCursorPosCallback(window, glfwCursorPosCallback);
        glfwSetMouseButtonCallback(window, glfwMouseButtonCallback);

        // 1. 初始化基础上下文
        VulkanContext vulkanContext;
        uint32_t glfwExtensionCount = 0;
        const char** glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);
        std::vector<const char*> requiredExtensions(glfwExtensions, glfwExtensions + glfwExtensionCount);
        vulkanContext.Initialize("Lizeral Sandbox", requiredExtensions);
        instance = (VkInstance)vulkanContext.GetNativeInstance();

        if (glfwCreateWindowSurface(instance, window, nullptr, &surface) != VK_SUCCESS) 
            throw std::runtime_error("Failed to create window surface!");

        {
            VulkanDevice vulkanDevice(&vulkanContext, surface);
            VulkanRenderer renderer(&vulkanContext, &vulkanDevice, window);

            // =======================================================
            // 2. ECS 装配阶段 (保持原样)
            // =======================================================
            Registry registry;
            CameraControlSystem cameraControlSystem;
            CameraSystem cameraSystem;

            Entity cameraEntity = registry.create();
            {
                auto& t = registry.emplace<TransformComponent>(cameraEntity);
                t.setPosition(Vector3(0.0f, 2.0f, 2.0f)); 
                auto& cam = registry.emplace<CameraComponent>(cameraEntity);
                cam.setFov(45.0f); cam.setAspect((float)WIDTH / (float)HEIGHT);
                cam.setzNear(0.001f); cam.setzFar(1000.0f); cam.setMain(true);
                auto& ctrl = registry.emplace<CameraControlComponent>(cameraEntity);
                ctrl.setMoveSpeed(3.0f); ctrl.setSensitivityX(0.1f); ctrl.setSensitivityY(0.1f);
                ctrl.setSpeedMutiplier(10.0f);
            }

            auto roomEntity = registry.create();
            auto& roomTrans = registry.emplace<TransformComponent>(roomEntity);
            roomTrans.setPosition(Vector3(0.0f, -1.0f, 10.0f));
            roomTrans.setScale(Vector3(0.05f, 0.05f, 0.05f)); 
            registry.emplace<VulkanModelComponent>(roomEntity).setVulkanModelPath("C:/Lizeral Engine/LizeralEngine0.0.1/asset/scene_without_window.glb"); 

            auto sponza = registry.create();
            auto& sponzaTrans = registry.emplace<TransformComponent>(sponza);
            sponzaTrans.setPosition(Vector3(100.0f, -1.0f, 10.0f));
            sponzaTrans.setScale(Vector3(1.0f, 1.1f, 1.1f)); 
            registry.emplace<VulkanModelComponent>(sponza).setVulkanModelPath("C:/Lizeral Engine/LizeralEngine0.0.1/asset/Sponza/Sponza/glTF/Sponza.gltf"); 

            // =======================================================
            // 3. ★ 核心：初始化新的渲染系统 ★
            // =======================================================
            VulkanRenderingSystem renderSystem;
            renderSystem.Initialize(&vulkanContext, &vulkanDevice, &renderer, WIDTH, HEIGHT);
            renderSystem.SetViewport(0, 0, WIDTH, HEIGHT);

            std::cout << "\n[Sandbox] Dynamic Rendering & Ray Tracing Ready! Hold RMB + WASD to move." << std::endl;
            float lastTime = glfwGetTime();

            // ====================================================================
            // 4. ★ 主循环：极致简化 ★
            // ====================================================================
            while (!glfwWindowShouldClose(window)) {
                
                float currentTime = glfwGetTime(); 
                float dt = currentTime - lastTime; 
                lastTime = currentTime;

                // 处理输入与物理/逻辑系统
                Input::GetInstance().Tick(); 
                if (Input::GetInstance().GetKey(Key::ESC)) glfwSetWindowShouldClose(window, true);
                glfwPollEvents();

                cameraControlSystem.Tick(dt, registry); 
                cameraSystem.Tick(registry); 

                // --- 将所有 Vulkan 渲染的脏活累活丢给系统 ---
                renderSystem.Tick(registry, dt);
            }

            std::cout << "--- EXITING PROGRAM ---" << std::endl;
            vkDeviceWaitIdle(vulkanDevice.GetNativeDevice());
            
            // 手动调用 Shutdown 确保在 m_device 被销毁前清理完毕
            renderSystem.Shutdown();
        } 
        
        vkDestroySurfaceKHR(instance, surface, nullptr);
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL ERROR] " << e.what() << std::endl;
    }

    if (window) {
        glfwDestroyWindow(window);
        glfwTerminate();
    }
    return 0;
}