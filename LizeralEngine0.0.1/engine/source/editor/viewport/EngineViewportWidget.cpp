#include "EngineViewportWidget.h"
#include <QApplication>
#include <QWindow>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#include <vulkan/vulkan_win32.h>
#endif

namespace Lizeral {

    EngineViewportWidget::EngineViewportWidget(Registry* registry, VulkanRenderingSystem* renderSys, QWidget* parent)
        : QWidget(parent), m_Registry(registry), m_RenderSystem(renderSys)
    {
        // 强制 Qt 交出绘制权
        setAttribute(Qt::WA_PaintOnScreen, true);
        setAttribute(Qt::WA_NativeWindow, true);
        setAttribute(Qt::WA_NoSystemBackground, true);
        setAttribute(Qt::WA_OpaquePaintEvent, true);
        setAutoFillBackground(false);

        setFocusPolicy(Qt::StrongFocus);
        setMouseTracking(true); // 开启鼠标追踪，这是漫游系统必需的

        m_resizeTimer = new QTimer(this);
        m_resizeTimer->setSingleShot(true); // 只触发一次
        connect(m_resizeTimer, &QTimer::timeout, this, &EngineViewportWidget::PerformResize);
    }

    EngineViewportWidget::~EngineViewportWidget() {
        CleanupVulkan();
    }

    // =========================================================================
    // Vulkan 生命周期与大小调整
    // =========================================================================
    void EngineViewportWidget::showEvent(QShowEvent* event) {
        QWidget::showEvent(event);
        if (!m_vulkanInitialized) {
            InitVulkan();
            m_vulkanInitialized = true;
        }
    }

    void EngineViewportWidget::InitVulkan() {
        m_context = std::make_unique<VulkanContext>();
        std::vector<const char*> windowExtensions = { "VK_KHR_surface", "VK_KHR_win32_surface" };
        m_context->Initialize("Lizeral Editor", windowExtensions);

#ifdef _WIN32
        HWND hwnd = reinterpret_cast<HWND>(winId());
        HINSTANCE hinstance = GetModuleHandle(nullptr);
        VkWin32SurfaceCreateInfoKHR createInfo{VK_STRUCTURE_TYPE_WIN32_SURFACE_CREATE_INFO_KHR};
        createInfo.hwnd = hwnd;
        createInfo.hinstance = hinstance;

        VkInstance vkInstance = static_cast<VkInstance>(m_context->GetNativeInstance());
        if (vkCreateWin32SurfaceKHR(vkInstance, &createInfo, nullptr, &m_surface) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create window surface from Qt Widget!");
        }
#endif

        m_device = std::make_unique<VulkanDevice>(m_context.get(), m_surface);
        m_renderer = std::make_unique<VulkanRenderer>(m_context.get(), m_device.get(), nullptr);
        
        // RenderSystem 由外部传入，我们只需调用其 Init
        if (m_RenderSystem) {
            m_RenderSystem->Initialize(m_context.get(), m_device.get(), m_renderer.get(), width(), height());
        }

        if (onInitVulkan) onInitVulkan();
    }

    void EngineViewportWidget::CleanupVulkan() {
        if (m_vulkanInitialized) {
            if (m_device) m_device->WaitIdle(); // 先死等所有任务做完
            
            // ==========================================
            // ★ 核心修复：必须先销毁 Renderer！
            // 这样它内部的 CommandBuffer 会被率先释放，
            // 从而解除对光追树 (TLAS) 的引用占用！
            // ==========================================
            m_renderer.reset(); 
            
            // 此时没有任何指令缓冲还在引用资源了，可以安全无痛地销毁渲染系统
            if (m_RenderSystem) m_RenderSystem->Shutdown();
            
            m_device.reset();   

            if (m_surface != VK_NULL_HANDLE && m_context) {
                VkInstance vkInstance = static_cast<VkInstance>(m_context->GetNativeInstance());
                vkDestroySurfaceKHR(vkInstance, m_surface, nullptr);
                m_surface = VK_NULL_HANDLE;
            }

            m_context->Shutdown();
            m_context.reset();  
            m_vulkanInitialized = false;
        }
    }

    void EngineViewportWidget::resizeEvent(QResizeEvent* event) {
        QWidget::resizeEvent(event);
        if (!m_vulkanInitialized || event->size().width() == 0) return;

        m_pendingSize = event->size();
        
        // 如果 200 毫秒内又触发了 resize，定时器会被重置
        // 只有当用户停止拖动 200 毫秒后，才会真正触发 PerformResize
        m_resizeTimer->start(200); 
    }

    // [+] 新增真正的重建函数
    void EngineViewportWidget::PerformResize() {
        m_device->WaitIdle(); 
        
        float dpiScale = this->devicePixelRatioF();
        int w = static_cast<int>(m_pendingSize.width() * dpiScale);
        int h = static_cast<int>(m_pendingSize.height() * dpiScale);
        if (w == 0 || h == 0) return;
        
        m_renderer->RecreateSwapchain(w, h);
        VkExtent2D actualExt = m_renderer->GetSwapchainExtent();
        
        if (m_RenderSystem) {
            m_RenderSystem->Shutdown();
            m_RenderSystem->Initialize(m_context.get(), m_device.get(), m_renderer.get(), actualExt.width, actualExt.height);
        }
    }

    // =========================================================================
    // 每帧渲染逻辑与“沉浸式”视口偏移
    // =========================================================================
    void EngineViewportWidget::paintEvent(QPaintEvent* event) {
        Q_UNUSED(event);
        if (!m_vulkanInitialized || !m_Registry || !m_RenderSystem) return;

        float dpiScale = this->devicePixelRatioF();
        
        // ★ 1. 获取最顶层的主窗口
        QWidget* mainWindow = this->window(); 
        
        int fullW = static_cast<int>(mainWindow->width() * dpiScale);
        int fullH = static_cast<int>(mainWindow->height() * dpiScale);
        if (fullH == 0) fullH = 1;

        // 2. 强制相机使用主窗口的长宽比，防止 UI 挤压变形
        float aspect = static_cast<float>(fullW) / static_cast<float>(fullH);
        auto view = m_Registry->view<CameraComponent>();
        for (auto entity : view) {
            auto& cam = view.get<CameraComponent>(entity);
            if (cam.isMain()) { 
                cam.setAspect(aspect);
                break;
            }
        }

        // ★ 3. 计算相对主窗口的物理正数偏移量
        QPoint topLeftPos = this->mapTo(mainWindow, QPoint(0, 0));
        int physicalOffsetX = static_cast<int>(topLeftPos.x() * dpiScale);
        int physicalOffsetY = static_cast<int>(topLeftPos.y() * dpiScale);

        // ★ 4. 将绝对真实偏移量传给 RenderingSystem
        m_RenderSystem->SetViewport(physicalOffsetX, physicalOffsetY, fullW, fullH);

        // 5. 驱动渲染
        m_RenderSystem->Tick(*m_Registry, 0.016f); 

        // 请求下一帧重绘，保持 60 帧率
        // update();
    }

    // =========================================================================
    // 完美的无限相机漫游系统
    // =========================================================================
    void EngineViewportWidget::mousePressEvent(QMouseEvent *event){
        if (event->button() == Qt::RightButton) {
            m_isRoaming = true;
            m_roamStartGlobalPos = QCursor::pos(); 
            m_virtualMousePos = event->pos(); 

            Input::GetInstance().SetMouseButtonDown(MouseButton::Right, true);
            setCursor(Qt::BlankCursor);
        }
    }

    void EngineViewportWidget::mouseReleaseEvent(QMouseEvent *event){
        if (event->button() == Qt::RightButton) {
            m_isRoaming = false;
            Input::GetInstance().SetMouseButtonDown(MouseButton::Right, false);
            
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

            Input::GetInstance().SetMousePosition(m_virtualMousePos.x(), m_virtualMousePos.y());
            QCursor::setPos(m_roamStartGlobalPos); // 强制把鼠标拉回锚点
        } else {
            Input::GetInstance().SetMousePosition(event->pos().x(), event->pos().y());
        }
    }

    void EngineViewportWidget::keyPressEvent(QKeyEvent *event){
        if (event->key() == Qt::Key_W) Input::GetInstance().SetKeyDown(Key::W, true);
        if (event->key() == Qt::Key_S) Input::GetInstance().SetKeyDown(Key::S, true);
        if (event->key() == Qt::Key_A) Input::GetInstance().SetKeyDown(Key::A, true);
        if (event->key() == Qt::Key_D) Input::GetInstance().SetKeyDown(Key::D, true);
        if (event->key() == Qt::Key_Q) Input::GetInstance().SetKeyDown(Key::Q, true);
        if (event->key() == Qt::Key_E) Input::GetInstance().SetKeyDown(Key::E, true);
        if (event->key() == Qt::Key_Shift) Input::GetInstance().SetKeyDown(Key::LEFT_SHIFT, true);
    }

    void EngineViewportWidget::keyReleaseEvent(QKeyEvent *event) {
        if (event->key() == Qt::Key_W) Input::GetInstance().SetKeyDown(Key::W, false);
        if (event->key() == Qt::Key_S) Input::GetInstance().SetKeyDown(Key::S, false);
        if (event->key() == Qt::Key_A) Input::GetInstance().SetKeyDown(Key::A, false);
        if (event->key() == Qt::Key_D) Input::GetInstance().SetKeyDown(Key::D, false);
        if (event->key() == Qt::Key_Q) Input::GetInstance().SetKeyDown(Key::Q, false);
        if (event->key() == Qt::Key_E) Input::GetInstance().SetKeyDown(Key::E, false);
        if (event->key() == Qt::Key_Shift) Input::GetInstance().SetKeyDown(Key::LEFT_SHIFT, false);
    }

} // namespace Lizeral