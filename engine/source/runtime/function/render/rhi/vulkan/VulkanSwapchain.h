// source/runtime/function/render/rhi/vulkan/VulkanSwapchain.h
#pragma once
#include <vulkan/vulkan.h>
#include <vector>

namespace Lizeral {

    class VulkanContext;
    class VulkanDevice;

    struct SwapChainSupportDetails {
        VkSurfaceCapabilitiesKHR capabilities;
        std::vector<VkSurfaceFormatKHR> formats;
        std::vector<VkPresentModeKHR> presentModes;
    };

    class VulkanSwapchain {
    public:
        VulkanSwapchain(VulkanContext* context, VulkanDevice* device, VkSurfaceKHR surface, uint32_t width, uint32_t height);
        ~VulkanSwapchain();

        VulkanSwapchain(const VulkanSwapchain&) = delete;
        VulkanSwapchain& operator=(const VulkanSwapchain&) = delete;

        VkSwapchainKHR GetNativeSwapchain() const { return m_swapchain; }
        VkFormat GetImageFormat() const { return m_imageFormat; }
        VkExtent2D GetExtent() const { return m_extent; }
        
        const std::vector<VkImageView>& GetImageViews() const { return m_imageViews; }
        const std::vector<VkImage>& GetNativeImages() const { return m_images; }

        VkImageLayout GetImageLayout(uint32_t index) const { return m_imageLayouts[index]; }
        void SetImageLayout(uint32_t index, VkImageLayout layout) { m_imageLayouts[index] = layout; }

    private:
        VulkanContext* m_context { nullptr };
        VulkanDevice* m_device { nullptr };
        VkSurfaceKHR m_surface { VK_NULL_HANDLE };

        VkSwapchainKHR m_swapchain { VK_NULL_HANDLE };
        VkFormat m_imageFormat;
        VkExtent2D m_extent;

        std::vector<VkImage> m_images;        
        std::vector<VkImageView> m_imageViews; 
        std::vector<VkImageLayout> m_imageLayouts;

        SwapChainSupportDetails querySwapChainSupport();
        VkSurfaceFormatKHR chooseSwapSurfaceFormat(const std::vector<VkSurfaceFormatKHR>& availableFormats);
        VkPresentModeKHR chooseSwapPresentMode(const std::vector<VkPresentModeKHR>& availablePresentModes);
        VkExtent2D chooseSwapExtent(const VkSurfaceCapabilitiesKHR& capabilities, uint32_t width, uint32_t height);

        void createSwapchain(uint32_t width, uint32_t height);
        void createImageViews();
    };

} // namespace Lizeral