#pragma once
#include <vulkan/vulkan.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <memory>

namespace Lizeral {

    class VulkanContext;
    class VulkanDevice;
    class VulkanSwapchain;
    class VulkanRenderPass;
    class VulkanFramebuffer;
    class VulkanCommandPool;
    class VulkanCommandBuffer;

    class VulkanRenderer {
    public:
        // 构造函数接管底层核心设备
        VulkanRenderer(VulkanContext* context, VulkanDevice* device, GLFWwindow* window);
        ~VulkanRenderer();

        // 禁用拷贝（渲染器是全局唯一的单例级别对象）
        VulkanRenderer(const VulkanRenderer&) = delete;
        VulkanRenderer& operator=(const VulkanRenderer&) = delete;

        // ====================================================================
        // ★ 核心生命周期 API (让外界告别繁琐的同步锁)
        // ====================================================================
        
        // 开启新的一帧：底层自动等 Fence、向 Swapchain 申请图片、开始录制 CommandBuffer
        // 返回的 CommandBuffer 将交给上层的 RenderGraph / 系统去录制绘制指令
        VkCommandBuffer BeginFrame();

        // 结束一帧：底层自动结束 CommandBuffer、提交给 GPU 队列，并推送到显示器
        void EndFrame();

        // ====================================================================
        // ★ Pass 接口层 (为将来的多 Pass 架构做准备)
        // ====================================================================
        
        // 开启向屏幕输出的“最终通道”。
        // 未来这里还会增加 BeginShadowPass, BeginGBufferPass 等，或者由 RenderGraph 接管
        void BeginSwapChainRenderPass(VkCommandBuffer commandBuffer);
        void EndSwapChainRenderPass(VkCommandBuffer commandBuffer);

        // ====================================================================
        // Getters (供外界的 PipelineBuilder 或 Material 系统获取必要状态)
        // ====================================================================
        bool IsFrameInProgress() const { return m_isFrameStarted; }
        VkCommandBuffer GetCurrentCommandBuffer() const;
        uint32_t GetCurrentFrameIndex() const { return m_currentFrame; }
        
        VulkanSwapchain* GetSwapchain() const { return m_swapchain.get(); }
        VulkanRenderPass* GetSwapChainRenderPass() const { return m_defaultRenderPass.get(); }

    private:
        // 内部初始化与清理逻辑
        void CreateCommandBuffers();
        void CreateSyncObjects();
        void RecreateSwapchain(); // 为将来窗口拖拽改变大小预留
        void FreeCommandBuffers();

        VulkanContext* m_context{ nullptr };
        VulkanDevice* m_device{ nullptr };
        GLFWwindow* m_window{ nullptr };

        // ★ 资源管理：使用 unique_ptr 自动管理生命周期
        std::unique_ptr<VulkanSwapchain> m_swapchain;
        std::unique_ptr<VulkanRenderPass> m_defaultRenderPass;
        std::unique_ptr<VulkanFramebuffer> m_defaultFramebuffers;
        std::unique_ptr<VulkanCommandPool> m_commandPool;
        std::vector<std::unique_ptr<VulkanCommandBuffer>> m_commandBuffers;

        // ★ 多帧并发 (Frames in Flight) 状态机
        const int MAX_FRAMES_IN_FLIGHT = 2;
        uint32_t m_currentFrame = 0;    // 当前是第几套 CPU 指挥系统 (0 或 1)
        uint32_t m_imageIndex = 0;      // 拿到了 Swapchain 的哪张图
        bool m_isFrameStarted = false;  // 防止业务层连续调用 BeginFrame 的状态锁

        // 并发同步原语
        std::vector<VkSemaphore> m_imageAvailableSemaphores;
        std::vector<VkSemaphore> m_renderFinishedSemaphores;
        std::vector<VkFence> m_inFlightFences;
    };

} // namespace Lizeral