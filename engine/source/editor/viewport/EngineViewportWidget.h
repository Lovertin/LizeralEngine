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

class EngineViewportWidget : public QWidget {
public:
    std::function<void()> onInitVulkan;
    
    EngineViewportWidget(Lizeral::Registry* registry, Lizeral::VulkanRenderingSystem* renderSys, QWidget* parent = nullptr);
    ~EngineViewportWidget();

    QPaintEngine* paintEngine() const override { return nullptr; } 

    void SetPhysicsSystem(Lizeral::PhysicsSystem* physSys) { m_physicsSystem = physSys; }

    void SetDebugLines(const std::vector<DebugLineVertex>& lines) {
            m_debugLines = lines;
    }

protected:

    void showEvent(QShowEvent* event) override;
    void paintEvent(QPaintEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

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

    std::unique_ptr<VulkanContext> m_context;
    std::unique_ptr<VulkanDevice> m_device;
    std::unique_ptr<VulkanRenderer> m_renderer;
    VkSurfaceKHR m_surface { VK_NULL_HANDLE };

    std::vector<DebugLineVertex> m_debugLines;

    bool m_isVulkanInitialized { false };
    bool m_isRoaming { false };
    QPoint m_roamStartGlobalPos; 
    QPointF m_virtualMousePos;
};

} // namespace Lizeral