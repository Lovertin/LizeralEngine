// source/runtime/function/render/rhi/vulkan/VulkanSwapchain.h
#pragma once
#include <vulkan/vulkan.h>
#include <vector>

namespace Lizeral {

    class VulkanContext;
    class VulkanDevice;

    // 辅助结构体：用来存放我们“查户口”查出来的支持细节
    struct SwapChainSupportDetails {
        VkSurfaceCapabilitiesKHR capabilities;
        std::vector<VkSurfaceFormatKHR> formats;
        std::vector<VkPresentModeKHR> presentModes;
    };

    class VulkanSwapchain {
    public:
        // 构造函数需要 Context(拿物理设备查户口)，Device(负责实际创建)，Surface(目标表面)，以及窗口当前的宽高
        VulkanSwapchain(VulkanContext* context, VulkanDevice* device, VkSurfaceKHR surface, uint32_t width, uint32_t height);
        ~VulkanSwapchain();

        // 禁用拷贝
        VulkanSwapchain(const VulkanSwapchain&) = delete;
        VulkanSwapchain& operator=(const VulkanSwapchain&) = delete;

        // --- 核心获取接口 ---
        VkSwapchainKHR GetNativeSwapchain() const { return m_swapchain; }
        VkFormat GetImageFormat() const { return m_imageFormat; }
        VkExtent2D GetExtent() const { return m_extent; }
        
        // 获取所有轮播图的“相框”，之后我们在渲染时需要把它们绑定到 Framebuffer 上
        const std::vector<VkImageView>& GetImageViews() const { return m_imageViews; }
        const std::vector<VkImage>& GetNativeImages() const { return m_images; }

    private:
        VulkanContext* m_context { nullptr };
        VulkanDevice* m_device { nullptr };
        VkSurfaceKHR m_surface { VK_NULL_HANDLE };

        VkSwapchainKHR m_swapchain { VK_NULL_HANDLE };
        VkFormat m_imageFormat;
        VkExtent2D m_extent;

        std::vector<VkImage> m_images;         // 交换链里的原生图片
        std::vector<VkImageView> m_imageViews; // 图片的视图（相框）

        // 内部初始化流程
        SwapChainSupportDetails querySwapChainSupport();
        VkSurfaceFormatKHR chooseSwapSurfaceFormat(const std::vector<VkSurfaceFormatKHR>& availableFormats);
        VkPresentModeKHR chooseSwapPresentMode(const std::vector<VkPresentModeKHR>& availablePresentModes);
        VkExtent2D chooseSwapExtent(const VkSurfaceCapabilitiesKHR& capabilities, uint32_t width, uint32_t height);

        void createSwapchain(uint32_t width, uint32_t height);
        void createImageViews();
    };

} // namespace Lizeral