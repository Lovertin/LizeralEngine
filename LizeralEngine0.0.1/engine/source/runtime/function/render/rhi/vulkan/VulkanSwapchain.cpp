// source/runtime/function/render/rhi/vulkan/VulkanSwapchain.cpp
#include "VulkanSwapchain.h"
#include "VulkanContext.h"
#include "VulkanDevice.h"
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <iostream>

namespace Lizeral {

    VulkanSwapchain::VulkanSwapchain(VulkanContext* context, VulkanDevice* device, VkSurfaceKHR surface, uint32_t width, uint32_t height)
        : m_context(context), m_device(device), m_surface(surface) {
        
        createSwapchain(width, height);
        createImageViews();
    }

    VulkanSwapchain::~VulkanSwapchain() {
        VkDevice device = m_device->GetNativeDevice();

        // 先砸碎所有的相框 (ImageView)
        for (auto imageView : m_imageViews) {
            vkDestroyImageView(device, imageView, nullptr);
        }

        // 再销毁轮播图机器 (Swapchain)
        if (m_swapchain != VK_NULL_HANDLE) {
            vkDestroySwapchainKHR(device, m_swapchain, nullptr);
            std::cout << "[VulkanSwapchain] Swapchain destroyed." << std::endl;
        }
    }

    SwapChainSupportDetails VulkanSwapchain::querySwapChainSupport() {
        SwapChainSupportDetails details;
        VkPhysicalDevice device = m_context->GetPhysicalDevice();

        // 1. 查基础能力 (比如最小/最大分辨率)
        vkGetPhysicalDeviceSurfaceCapabilitiesKHR(device, m_surface, &details.capabilities);

        // 2. 查支持的色彩格式
        uint32_t formatCount;
        vkGetPhysicalDeviceSurfaceFormatsKHR(device, m_surface, &formatCount, nullptr);
        if (formatCount != 0) {
            details.formats.resize(formatCount);
            vkGetPhysicalDeviceSurfaceFormatsKHR(device, m_surface, &formatCount, details.formats.data());
        }

        // 3. 查支持的呈现模式
        uint32_t presentModeCount;
        vkGetPhysicalDeviceSurfacePresentModesKHR(device, m_surface, &presentModeCount, nullptr);
        if (presentModeCount != 0) {
            details.presentModes.resize(presentModeCount);
            vkGetPhysicalDeviceSurfacePresentModesKHR(device, m_surface, &presentModeCount, details.presentModes.data());
        }

        return details;
    }

    VkSurfaceFormatKHR VulkanSwapchain::chooseSwapSurfaceFormat(const std::vector<VkSurfaceFormatKHR>& availableFormats) {
        // 首选 SRGB 颜色空间（人眼看着最舒服，也是主流 PBR 材质的基础）
        for (const auto& availableFormat : availableFormats) {
            if (availableFormat.format == VK_FORMAT_B8G8R8A8_SRGB && availableFormat.colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR) {
                return availableFormat;
            }
        }
        // 如果没有首选，就随便拿第一个
        return availableFormats[0];
    }

    VkPresentModeKHR VulkanSwapchain::chooseSwapPresentMode(const std::vector<VkPresentModeKHR>& availablePresentModes) {
        // 尝试寻找 Mailbox（三重缓冲，延迟极低，画面不撕裂，适合游戏引擎）
        for (const auto& availablePresentMode : availablePresentModes) {
            if (availablePresentMode == VK_PRESENT_MODE_MAILBOX_KHR) {
                return availablePresentMode;
            }
        }
        // 如果没有，就退回到绝对都有的 FIFO (传统的垂直同步 V-Sync)
        return VK_PRESENT_MODE_FIFO_KHR;
    }

    VkExtent2D VulkanSwapchain::chooseSwapExtent(const VkSurfaceCapabilitiesKHR& capabilities, uint32_t width, uint32_t height) {
        if (capabilities.currentExtent.width != std::numeric_limits<uint32_t>::max()) {
            return capabilities.currentExtent;
        } else {
            // 如果窗口管理器允许我们自己定大小，我们就把传进来的宽高裁剪到显卡允许的范围内
            VkExtent2D actualExtent = { width, height };
            actualExtent.width = std::clamp(actualExtent.width, capabilities.minImageExtent.width, capabilities.maxImageExtent.width);
            actualExtent.height = std::clamp(actualExtent.height, capabilities.minImageExtent.height, capabilities.maxImageExtent.height);
            return actualExtent;
        }
    }

    void VulkanSwapchain::createSwapchain(uint32_t width, uint32_t height) {
        SwapChainSupportDetails swapChainSupport = querySwapChainSupport();

        VkSurfaceFormatKHR surfaceFormat = chooseSwapSurfaceFormat(swapChainSupport.formats);
        VkPresentModeKHR presentMode = chooseSwapPresentMode(swapChainSupport.presentModes);
        VkExtent2D extent = chooseSwapExtent(swapChainSupport.capabilities, width, height);

        // 决定轮播图里有几张图片（通常推荐 最小要求数 + 1，实现三重缓冲）
        uint32_t imageCount = swapChainSupport.capabilities.minImageCount + 1;
        if (swapChainSupport.capabilities.maxImageCount > 0 && imageCount > swapChainSupport.capabilities.maxImageCount) {
            imageCount = swapChainSupport.capabilities.maxImageCount;
        }

        VkSwapchainCreateInfoKHR createInfo{};
        createInfo.sType = VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR;
        createInfo.surface = m_surface;
        createInfo.minImageCount = imageCount;
        createInfo.imageFormat = surfaceFormat.format;
        createInfo.imageColorSpace = surfaceFormat.colorSpace;
        createInfo.imageExtent = extent;
        createInfo.imageArrayLayers = 1;
        createInfo.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT; // 我们要往这些图片里画颜色

        // 队列家族的处理：如果画图车间和呈现车间不一样，我们需要让图片可以在两个车间流通
        uint32_t indices[] = { m_device->GetGraphicsQueueFamily(), m_device->GetPresentQueueFamily()};
        if (indices[0] != indices[1]) {
            createInfo.imageSharingMode = VK_SHARING_MODE_CONCURRENT;
            createInfo.queueFamilyIndexCount = 2;
            createInfo.pQueueFamilyIndices = indices;
        } else {
            createInfo.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE; // 大部分情况走这里（性能更好）
        }

        createInfo.preTransform = swapChainSupport.capabilities.currentTransform;
        createInfo.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR; // 不透明窗口
        createInfo.presentMode = presentMode;
        createInfo.clipped = VK_TRUE; // 开启裁剪，被其他窗口挡住的像素就不画了，提升性能
        createInfo.oldSwapchain = VK_NULL_HANDLE;

        if (vkCreateSwapchainKHR(m_device->GetNativeDevice(), &createInfo, nullptr, &m_swapchain) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create swapchain!");
        }

        // 记录格式和分辨率，供后续渲染使用
        m_imageFormat = surfaceFormat.format;
        m_extent = extent;

        // 获取显卡自动帮我们分配好的那几张图片的把柄
        vkGetSwapchainImagesKHR(m_device->GetNativeDevice(), m_swapchain, &imageCount, nullptr);
        m_images.resize(imageCount);
        vkGetSwapchainImagesKHR(m_device->GetNativeDevice(), m_swapchain, &imageCount, m_images.data());

        std::cout << "[VulkanSwapchain] Swapchain created successfully! Images count: " << imageCount << std::endl;
    }

    void VulkanSwapchain::createImageViews() {
        m_imageViews.resize(m_images.size());

        for (size_t i = 0; i < m_images.size(); i++) {
            VkImageViewCreateInfo createInfo{};
            createInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
            createInfo.image = m_images[i];
            createInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
            createInfo.format = m_imageFormat;

            // RGBA 通道映射（保持原样）
            createInfo.components.r = VK_COMPONENT_SWIZZLE_IDENTITY;
            createInfo.components.g = VK_COMPONENT_SWIZZLE_IDENTITY;
            createInfo.components.b = VK_COMPONENT_SWIZZLE_IDENTITY;
            createInfo.components.a = VK_COMPONENT_SWIZZLE_IDENTITY;

            // 告诉 Vulkan 这是一张纯颜色贴图，没有 Mipmap，只有 1 层
            createInfo.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
            createInfo.subresourceRange.baseMipLevel = 0;
            createInfo.subresourceRange.levelCount = 1;
            createInfo.subresourceRange.baseArrayLayer = 0;
            createInfo.subresourceRange.layerCount = 1;

            if (vkCreateImageView(m_device->GetNativeDevice(), &createInfo, nullptr, &m_imageViews[i]) != VK_SUCCESS) {
                throw std::runtime_error("Failed to create image views!");
            }
        }
        std::cout << "[VulkanSwapchain] Created " << m_imageViews.size() << " Image Views." << std::endl;
    }

} // namespace Lizeral