#pragma once

#include <vulkan/vulkan.h>
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

        VulkanRenderer(VulkanContext* context, VulkanDevice* device, uint32_t width,uint32_t height);
        ~VulkanRenderer();

        VkCommandBuffer BeginFrame();
        void EndFrame();

        void BeginRendering(VkCommandBuffer commandBuffer);
        void EndRendering(VkCommandBuffer commandBuffer);

        VkFormat GetSwapchainFormat() const;

        VkExtent2D GetSwapchainExtent() const; 

        void RecreateSwapchain(uint32_t width, uint32_t height);
        
        bool IsSwapchainOutdated() const { return m_swapchainOutdated; }

    private:
        void CreateCommandBuffers();
        void FreeCommandBuffers();
        void CreateSyncObjects();
        void CleanupDepthResources();

        VkCommandBuffer GetCurrentCommandBuffer() const;

        VulkanContext* m_context;
        VulkanDevice* m_device;
        uint32_t m_width = 0;
        uint32_t m_height = 0;
        bool m_swapchainOutdated = false;

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