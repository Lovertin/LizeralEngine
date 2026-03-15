#pragma once

#include <QWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <QTimer>
#include <memory>
#include <functional>

// ECS
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"

// Vulkan
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"
#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystem.h"

// 物理与输入
// #include "runtime/function/physics/PhysicsSystem.h" // 替换为你的真实路径
#include "runtime/function/input/input.h"

namespace Lizeral {

    class EngineViewportWidget : public QWidget {
        Q_OBJECT

    public:
        std::function<void()> onInitVulkan;

        explicit EngineViewportWidget(Registry* registry, VulkanRenderingSystem* renderSys, QWidget* parent = nullptr);
        ~EngineViewportWidget() override;

        // 设置物理系统
        // void SetPhysicsSystem(PhysicsSystem* physSys) { m_physicsSystem = physSys; }

    protected:
        void showEvent(QShowEvent* event) override;
        void resizeEvent(QResizeEvent* event) override;
        void paintEvent(QPaintEvent* event) override;

        void mousePressEvent(QMouseEvent *event) override;
        void mouseReleaseEvent(QMouseEvent *event) override;
        void mouseMoveEvent(QMouseEvent *event) override;
        void keyPressEvent(QKeyEvent *event) override;
        void keyReleaseEvent(QKeyEvent *event) override;
        QPaintEngine* paintEngine() const override { return nullptr; } 
        

    private slots:
        void PerformResize(); // [+] 新增槽函数

    private:
        void InitVulkan();
        void CleanupVulkan();

        bool m_vulkanInitialized = false;

        // --- Vulkan 核心组件 ---
        std::unique_ptr<VulkanContext> m_context;
        VkSurfaceKHR m_surface { VK_NULL_HANDLE };
        std::unique_ptr<VulkanDevice> m_device;
        std::unique_ptr<VulkanRenderer> m_renderer;
        
        // 引用外部传入的系统
        Registry* m_Registry { nullptr };
        VulkanRenderingSystem* m_RenderSystem { nullptr };
        // PhysicsSystem* m_physicsSystem { nullptr };

        // --- 漫游系统状态 ---
        bool m_isRoaming { false };
        QPoint m_roamStartGlobalPos; 
        QPointF m_virtualMousePos;

        QTimer* m_resizeTimer { nullptr }; // [+] 新增定时器
        QSize m_pendingSize;               // [+] 记录最新的窗口大小
    };

} // namespace Lizeral