#include "EngineViewportWidget.h"
#include <QWindow>
#include <QGuiApplication>

#include <windows.h>

#define VK_USE_PLATFORM_WIN32_KHR
#include <vulkan/vulkan_win32.h>

namespace Lizeral {

    EngineViewportWidget::EngineViewportWidget(Lizeral::Registry* registry, Lizeral::VulkanRenderingSystem* renderSys, QWidget* parent) 
        : QWidget(parent), m_Registry(registry), m_RenderSystem(renderSys) 
    {
        setFocusPolicy(Qt::StrongFocus);
        setMouseTracking(true); 
        setAttribute(Qt::WA_PaintOnScreen);
        setAttribute(Qt::WA_NoSystemBackground);
        setAttribute(Qt::WA_NativeWindow);
        setAttribute(Qt::WA_OpaquePaintEvent);
    }

    EngineViewportWidget::~EngineViewportWidget() {
        if (m_isVulkanInitialized) {
            // wait GPU idle
            vkDeviceWaitIdle(m_device->GetNativeDevice());
            
            m_renderer.reset();
            
            m_RenderSystem->Shutdown();
            
            m_device.reset();
            if (m_surface != VK_NULL_HANDLE) {
                vkDestroySurfaceKHR(static_cast<VkInstance>(m_context->GetNativeInstance()), m_surface, nullptr);
            }
            m_context.reset();
        }
    }

    void EngineViewportWidget::showEvent(QShowEvent* event) {
        QWidget::showEvent(event);
        //Initialize Vulkan
        if (!m_isVulkanInitialized) {
            InitVulkan();
        }
    }

    void EngineViewportWidget::InitVulkan() {
        m_context = std::make_unique<VulkanContext>();
        std::vector<const char*> extensions = { "VK_KHR_surface", "VK_KHR_win32_surface", "VK_EXT_debug_utils" };
        m_context->Initialize("Lizeral Editor", extensions);
        VkInstance instance = static_cast<VkInstance>(m_context->GetNativeInstance());

        VkWin32SurfaceCreateInfoKHR createInfo{VK_STRUCTURE_TYPE_WIN32_SURFACE_CREATE_INFO_KHR};
        createInfo.hwnd = reinterpret_cast<HWND>(this->winId());
        createInfo.hinstance = GetModuleHandle(nullptr);
        if (vkCreateWin32SurfaceKHR(instance, &createInfo, nullptr, &m_surface) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create Win32 Vulkan Surface!");
        }

        m_device = std::make_unique<VulkanDevice>(m_context.get(), m_surface);

        float dpiScale = this->devicePixelRatioF();
        uint32_t w = static_cast<uint32_t>(this->width() * dpiScale);
        uint32_t h = static_cast<uint32_t>(this->height() * dpiScale);

        m_renderer = std::make_unique<VulkanRenderer>(m_context.get(), m_device.get(), w, h);
        VkExtent2D actualExt = m_renderer->GetSwapchainExtent();
        m_RenderSystem->Initialize(m_context.get(), m_device.get(), m_renderer.get(), actualExt.width, actualExt.height);
        
        m_isVulkanInitialized = true;

        if (onInitVulkan) {
            onInitVulkan(); 
        }

        this->update();
    }

    void EngineViewportWidget::resizeEvent(QResizeEvent* event) {
        QWidget::resizeEvent(event);
        if (!m_isVulkanInitialized || event->size().width() == 0 || event->size().height() == 0) return;

        float dpiScale = this->devicePixelRatioF();
        uint32_t w = static_cast<uint32_t>(event->size().width() * dpiScale);
        uint32_t h = static_cast<uint32_t>(event->size().height() * dpiScale);
        
        if (w > 0 && h > 0) {
            VkExtent2D currentExt = m_renderer->GetSwapchainExtent();
            if (currentExt.width == w && currentExt.height == h && !m_renderer->IsSwapchainOutdated()) {
                return;
            }
            m_renderer->RecreateSwapchain(w, h);
            VkExtent2D actualExt = m_renderer->GetSwapchainExtent();
            m_RenderSystem->Resize(actualExt.width, actualExt.height);
        }
    }

    void EngineViewportWidget::paintEvent(QPaintEvent* event) {
        if (!m_isVulkanInitialized || !m_Registry) return;

        if (m_renderer->IsSwapchainOutdated()) {
            float dpiScale = this->devicePixelRatioF();
            uint32_t w = static_cast<uint32_t>(this->width() * dpiScale);
            uint32_t h = static_cast<uint32_t>(this->height() * dpiScale);
            if (w > 0 && h > 0) {
                m_renderer->RecreateSwapchain(w, h);
                VkExtent2D actualExt = m_renderer->GetSwapchainExtent();
                m_RenderSystem->Resize(actualExt.width, actualExt.height);
            }
        }

        VkExtent2D currentExt = m_renderer->GetSwapchainExtent();
        m_RenderSystem->SetViewport(0, 0, currentExt.width, currentExt.height);

        // Invoke the black box mechanism
        m_RenderSystem->Tick(*m_Registry, 0.016f,m_debugLines); 
    }

    void EngineViewportWidget::mousePressEvent(QMouseEvent *event){
        setFocus(Qt::MouseFocusReason);

        if (event->button() == Qt::RightButton) {
            m_isRoaming = true;
            m_roamStartGlobalPos = QCursor::pos(); 
            m_virtualMousePos = event->pos(); 
            Lizeral::Input::GetInstance().SetMouseButtonDown(Lizeral::MouseButton::Right, true);
            Lizeral::Input::GetInstance().ResetMouse();
            setCursor(Qt::BlankCursor);
        } else if (event->button() == Qt::LeftButton) {
            // Reset first-sample mouse delta when returning focus to viewport.
            Lizeral::Input::GetInstance().ResetMouse();
            Lizeral::Input::GetInstance().SetMousePosition(event->pos().x(), event->pos().y());
        }
    }

    void EngineViewportWidget::mouseReleaseEvent(QMouseEvent *event){
        if (event->button() == Qt::RightButton) {
            m_isRoaming = false;
            Lizeral::Input::GetInstance().SetMouseButtonDown(Lizeral::MouseButton::Right, false);
            setCursor(Qt::ArrowCursor); 
            QCursor::setPos(m_roamStartGlobalPos);
        }
    }

    void EngineViewportWidget::mouseMoveEvent(QMouseEvent *event){
        if (m_isRoaming) {
            QPoint currentGlobalPos = QCursor::pos();
            if (currentGlobalPos == m_roamStartGlobalPos) return;

            int dx = currentGlobalPos.x() - m_roamStartGlobalPos.x();
            int dy = currentGlobalPos.y() - m_roamStartGlobalPos.y();

            m_virtualMousePos.setX(m_virtualMousePos.x() + dx);
            m_virtualMousePos.setY(m_virtualMousePos.y() + dy);

            Lizeral::Input::GetInstance().SetMousePosition(m_virtualMousePos.x(), m_virtualMousePos.y());
            QCursor::setPos(m_roamStartGlobalPos);
        } else {
            Lizeral::Input::GetInstance().SetMousePosition(event->pos().x(), event->pos().y());
        }
    }

    void EngineViewportWidget::keyPressEvent(QKeyEvent *event){
        if (event->key() == Qt::Key_W) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::W, true);
        if (event->key() == Qt::Key_S) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::S, true);
        if (event->key() == Qt::Key_A) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::A, true);
        if (event->key() == Qt::Key_D) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::D, true);
        if (event->key() == Qt::Key_Q) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::Q, true);
        if (event->key() == Qt::Key_E) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::E, true);
    }

    void EngineViewportWidget::keyReleaseEvent(QKeyEvent *event) {
        if (event->key() == Qt::Key_W) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::W, false);
        if (event->key() == Qt::Key_S) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::S, false);
        if (event->key() == Qt::Key_A) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::A, false);
        if (event->key() == Qt::Key_D) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::D, false);
        if (event->key() == Qt::Key_Q) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::Q, false);
        if (event->key() == Qt::Key_E) Lizeral::Input::GetInstance().SetKeyDown(Lizeral::Key::E, false);
    }

    void EngineViewportWidget::focusInEvent(QFocusEvent* event) {
        QWidget::focusInEvent(event);
        // Avoid a large first mouse delta when focus returns to viewport.
        Lizeral::Input::GetInstance().ResetMouse();
    }

    void EngineViewportWidget::focusOutEvent(QFocusEvent* event) {
        QWidget::focusOutEvent(event);

        // Leaving viewport focus (e.g. clicking editor panels) should not keep stale input states.
        m_isRoaming = false;
        setCursor(Qt::ArrowCursor);
        Lizeral::Input::GetInstance().ResetState();
    }
}
