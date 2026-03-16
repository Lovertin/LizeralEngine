#pragma once
#include "editor/EditorHeader.h"
#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystem.h"
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"

#include <QWidget>
#include <QEvent>
#include <memory>
#include <vulkan/vulkan.h>

namespace Lizeral {

// ★ 核心改变：不再继承 QOpenGLWidget，而是直接继承原生的 QWidget
class EngineViewportWidget : public QWidget {
public:
    std::function<void()> onInitVulkan;
    
    EngineViewportWidget(Lizeral::Registry* registry, Lizeral::VulkanRenderingSystem* renderSys, QWidget* parent = nullptr);
    ~EngineViewportWidget();

    // ★ 核心魔法：返回 nullptr 彻底禁用 Qt 自己的绘制引擎，让 Vulkan 独占这块屏幕像素
    QPaintEngine* paintEngine() const override { return nullptr; } 

    void SetPhysicsSystem(Lizeral::PhysicsSystem* physSys) { m_physicsSystem = physSys; }

protected:
    // 生命周期与渲染事件
    void showEvent(QShowEvent* event) override;
    void paintEvent(QPaintEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

    // 键鼠输入事件 (保持不变)
    void mousePressEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;

private:
    void InitVulkan();

    Lizeral::Registry* m_Registry { nullptr };
    Lizeral::VulkanRenderingSystem* m_RenderSystem { nullptr };
    Lizeral::PhysicsSystem* m_physicsSystem {nullptr};

    // --- Vulkan 核心 RHI 实例 (由 Widget 持有生命周期) ---
    std::unique_ptr<VulkanContext> m_context;
    std::unique_ptr<VulkanDevice> m_device;
    std::unique_ptr<VulkanRenderer> m_renderer;
    VkSurfaceKHR m_surface { VK_NULL_HANDLE };

    bool m_isVulkanInitialized { false };
    bool m_isRoaming { false };
    QPoint m_roamStartGlobalPos; 
    QPointF m_virtualMousePos;
};

} // namespace Lizeral