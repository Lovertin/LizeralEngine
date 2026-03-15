#pragma once

#include <vulkan/vulkan.h>
#include <GLFW/glfw3.h>
#include <memory>
#include <vector>

namespace Lizeral {

    class VulkanContext;
    class VulkanDevice;
    class VulkanSwapchain;
    class VulkanCommandPool;
    class VulkanCommandBuffer;

    class VulkanRenderer {
    public:
        VulkanRenderer(VulkanContext* context, VulkanDevice* device, GLFWwindow* window);
        ~VulkanRenderer();

        VulkanRenderer(const VulkanRenderer&) = delete;
        VulkanRenderer& operator=(const VulkanRenderer&) = delete;

        // --- 核心渲染循环接口 ---
        VkCommandBuffer BeginFrame();
        void EndFrame();

        // ★ 动态渲染核心：随叫随到，指哪画哪！
        void BeginRendering(VkCommandBuffer commandBuffer);
        void EndRendering(VkCommandBuffer commandBuffer);

        // --- 获取状态接口 ---
        VkCommandBuffer GetCurrentCommandBuffer() const;
        int GetFrameIndex() const { return m_currentFrame; }
        bool IsFrameInProgress() const { return m_isFrameStarted; }
        
        // 供管线创建时获取颜色格式
        VkFormat GetSwapchainFormat() const; 

        VkExtent2D GetSwapchainExtent() const ;

        void RecreateSwapchain(int width=1280, int height=720);

    private:
        VulkanContext* m_context { nullptr };
        VulkanDevice* m_device { nullptr };
        GLFWwindow* m_window { nullptr };

        std::unique_ptr<VulkanSwapchain> m_swapchain;
        std::unique_ptr<VulkanCommandPool> m_commandPool;
        std::vector<std::unique_ptr<VulkanCommandBuffer>> m_commandBuffers;

        // ★ 深度图资源 (初始化为空)
        VkImage m_depthImage { VK_NULL_HANDLE };
        VkDeviceMemory m_depthImageMemory { VK_NULL_HANDLE };
        VkImageView m_depthImageView { VK_NULL_HANDLE };

        // 同步对象
        std::vector<VkSemaphore> m_imageAvailableSemaphores;
        std::vector<VkSemaphore> m_renderFinishedSemaphores;
        std::vector<VkFence> m_inFlightFences;

        uint32_t m_imageIndex { 0 };
        int m_currentFrame { 0 };
        bool m_isFrameStarted { false };
        
        const int MAX_FRAMES_IN_FLIGHT = 2;

        // 内部流程函数
        void CreateCommandBuffers();
        void FreeCommandBuffers();
        void CreateSyncObjects();
        void CleanupDepthResources();
    };

} // namespace Lizeral