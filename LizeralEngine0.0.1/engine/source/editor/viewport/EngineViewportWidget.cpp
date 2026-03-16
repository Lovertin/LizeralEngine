#include "EngineViewportWidget.h"
#include <QWindow>
#include <QGuiApplication>

// ★ 核心修复：必须先包含 Windows 头文件，再包含 Vulkan Win32 头文件
// 加入这两个宏可以防止 Windows API 污染 C++ 标准库（比如 min/max 冲突）
#include <windows.h>

#define VK_USE_PLATFORM_WIN32_KHR
#include <vulkan/vulkan_win32.h>

namespace Lizeral {

    EngineViewportWidget::EngineViewportWidget(Lizeral::Registry* registry, Lizeral::VulkanRenderingSystem* renderSys, QWidget* parent) 
        : QWidget(parent), m_Registry(registry), m_RenderSystem(renderSys) 
    {
        setFocusPolicy(Qt::StrongFocus);
        setMouseTracking(true); 

        // 告诉操作系统：我们自己接管这个窗口的绘制，不要用系统的默认背景去刷它（防止闪烁）
        setAttribute(Qt::WA_PaintOnScreen);
        setAttribute(Qt::WA_NoSystemBackground);
        setAttribute(Qt::WA_NativeWindow);
        setAttribute(Qt::WA_OpaquePaintEvent);
    }

    EngineViewportWidget::~EngineViewportWidget() {
        if (m_isVulkanInitialized) {
            // 等待 GPU 停下
            vkDeviceWaitIdle(m_device->GetNativeDevice());
            
            m_renderer.reset();
            
            // 现在可以极其安全地销毁 RenderingSystem 里的光追树、管线和图片了
            m_RenderSystem->Shutdown();
            
            // 最后释放设备和上下文
            m_device.reset();
            if (m_surface != VK_NULL_HANDLE) {
                vkDestroySurfaceKHR(static_cast<VkInstance>(m_context->GetNativeInstance()), m_surface, nullptr);
            }
            m_context.reset();
        }
    }

    void EngineViewportWidget::showEvent(QShowEvent* event) {
        QWidget::showEvent(event);
        // 当 UI 面板首次显示在屏幕上时，初始化 Vulkan
        if (!m_isVulkanInitialized) {
            InitVulkan();
        }
    }

    void EngineViewportWidget::InitVulkan() {
        // 1. 初始化 Context
        m_context = std::make_unique<VulkanContext>();
        std::vector<const char*> extensions = { "VK_KHR_surface", "VK_KHR_win32_surface", "VK_EXT_debug_utils" };
        m_context->Initialize("Lizeral Editor", extensions);
        VkInstance instance = static_cast<VkInstance>(m_context->GetNativeInstance());

        // 2. 利用 Qt 获取操作系统的 HWND (句柄)，强行绑定 Vulkan Surface！
        VkWin32SurfaceCreateInfoKHR createInfo{VK_STRUCTURE_TYPE_WIN32_SURFACE_CREATE_INFO_KHR};
        createInfo.hwnd = reinterpret_cast<HWND>(this->winId());
        createInfo.hinstance = GetModuleHandle(nullptr);
        if (vkCreateWin32SurfaceKHR(instance, &createInfo, nullptr, &m_surface) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create Win32 Vulkan Surface!");
        }

        // 3. 创建 Device
        m_device = std::make_unique<VulkanDevice>(m_context.get(), m_surface);

        // 4. 获取真实的物理 DPI 像素
        float dpiScale = this->devicePixelRatioF();
        uint32_t w = static_cast<uint32_t>(this->width() * dpiScale);
        uint32_t h = static_cast<uint32_t>(this->height() * dpiScale);

        // 5. 初始化解耦后的 Renderer 与 RenderingSystem
        m_renderer = std::make_unique<VulkanRenderer>(m_context.get(), m_device.get(), w, h);
        VkExtent2D actualExt = m_renderer->GetSwapchainExtent();
        m_RenderSystem->Initialize(m_context.get(), m_device.get(), m_renderer.get(), w, h);
        
        m_isVulkanInitialized = true;

        // 触发 Editor 加载资源 (回调到 LizeralEditorWindow.cpp 里的 populateTestData)
        if (onInitVulkan) {
            onInitVulkan(); 
        }

        this->update();
    }

    void EngineViewportWidget::resizeEvent(QResizeEvent* event) {
        QWidget::resizeEvent(event);
        if (!m_isVulkanInitialized || event->size().width() == 0 || event->size().height() == 0) return;

        // ★ 核心桥梁：拦截 Qt 的 Resize 事件，翻译成真正的物理像素塞给 RHI 黑盒！
        float dpiScale = this->devicePixelRatioF();
        uint32_t w = static_cast<uint32_t>(event->size().width() * dpiScale);
        uint32_t h = static_cast<uint32_t>(event->size().height() * dpiScale);
        
        if (w > 0 && h > 0) {
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

        if (QWidget* mainWindow = this->topLevelWidget()) {
            float dpiScale = this->devicePixelRatioF();
            // 获取外面那个最大的 Editor MainWindow 的尺寸
            uint32_t fullW = static_cast<uint32_t>(mainWindow->width() * dpiScale);
            uint32_t fullH = static_cast<uint32_t>(mainWindow->height() * dpiScale);
            
            // 把这套“全局视野尺寸”喂给系统，你的 Tick 里的投影矩阵就会算出完美的无畸变画面！
            m_RenderSystem->SetViewport(0, 0, fullW, fullH);
        }

        // 触发黑盒渲染
        m_RenderSystem->Tick(*m_Registry, 0.016f); 
    }

    // =========================================================================
    // 键鼠事件 (全部保持原样，直接透传给 Input 系统)
    // =========================================================================
    void EngineViewportWidget::mousePressEvent(QMouseEvent *event){
        if (event->button() == Qt::RightButton) {
            m_isRoaming = true;
            m_roamStartGlobalPos = QCursor::pos(); 
            m_virtualMousePos = event->pos(); 
            Lizeral::Input::GetInstance().SetMouseButtonDown(Lizeral::MouseButton::Right, true);
            setCursor(Qt::BlankCursor);
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
}