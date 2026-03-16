#pragma once

#include <vulkan/vulkan.h>
#include <GLFW/glfw3.h> // 恢复对 GLFW 的依赖
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
        static constexpr int MAX_FRAMES_IN_FLIGHT = 2;

        VulkanRenderer(VulkanContext* context, VulkanDevice* device, GLFWwindow* window);
        ~VulkanRenderer();

        VkCommandBuffer BeginFrame();
        void EndFrame();

        void BeginRendering(VkCommandBuffer commandBuffer);
        void EndRendering(VkCommandBuffer commandBuffer);

        VkFormat GetSwapchainFormat() const;
        
        // 获取真实的物理像素分辨率
        VkExtent2D GetSwapchainExtent() const; 

        // 恢复为原版：自动通过 glfw 获取真实分辨率
        void RecreateSwapchain();

    private:
        void CreateCommandBuffers();
        void FreeCommandBuffers();
        void CreateSyncObjects();
        void CleanupDepthResources();

        VkCommandBuffer GetCurrentCommandBuffer() const;

        VulkanContext* m_context;
        VulkanDevice* m_device;
        GLFWwindow* m_window; // 恢复 GLFWwindow

        std::unique_ptr<VulkanSwapchain> m_swapchain;
        std::unique_ptr<VulkanCommandPool> m_commandPool;
        std::vector<std::unique_ptr<VulkanCommandBuffer>> m_commandBuffers;

        VkImage m_depthImage { VK_NULL_HANDLE };
        VkDeviceMemory m_depthImageMemory { VK_NULL_HANDLE };
        VkImageView m_depthImageView { VK_NULL_HANDLE };

        std::vector<VkSemaphore> m_imageAvailableSemaphores;
        std::vector<VkSemaphore> m_renderFinishedSemaphores;
        std::vector<VkFence> m_inFlightFences;

        uint32_t m_currentFrame = 0;
        uint32_t m_imageIndex;
        bool m_isFrameStarted = false;
    };

} // namespace Lizeral