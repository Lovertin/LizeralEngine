// source/runtime/function/render/rhi/vulkan/VulkanDevice.h
#pragma once
#include "runtime/function/render/rhi/RHIDevice.h"
#include <vulkan/vulkan.h>
#include <vma/vk_mem_alloc.h> 
#include <optional>

namespace Lizeral {
    
    class VulkanContext;

    struct QueueFamilyIndices {
        std::optional<uint32_t> graphicsFamily;
        std::optional<uint32_t> presentFamily;
        
        bool isComplete() {
            return graphicsFamily.has_value();
        }
    };

    class VulkanDevice : public IRHIDevice {
    public:

        VulkanDevice(VulkanContext* context,VkSurfaceKHR surface);
        ~VulkanDevice() override;

        void WaitIdle() override;

        VkDevice GetNativeDevice() const { return m_logicalDevice; }
        VkQueue GetGraphicsQueue() const { return m_graphicsQueue; }
        VkQueue GetPresentQueue() const { return m_presentQueue; }
        VmaAllocator GetAllocator() const { return m_allocator; }
        uint32_t GetGraphicsQueueFamily() const { return m_queueIndices.graphicsFamily.value(); }
        uint32_t GetPresentQueueFamily() const {return m_queueIndices.presentFamily.value();}

        VulkanContext* GetContext() const { return m_context;}
        VkSurfaceKHR GetSurface() const { return m_surface; }

    private:
        VulkanContext* m_context { nullptr };
        VkSurfaceKHR m_surface{ VK_NULL_HANDLE };
        
        VkDevice m_logicalDevice { VK_NULL_HANDLE };
        VkQueue m_graphicsQueue { VK_NULL_HANDLE };
        VkQueue m_presentQueue { VK_NULL_HANDLE };
        VmaAllocator m_allocator { VK_NULL_HANDLE };
        
        QueueFamilyIndices m_queueIndices;

        void findQueueFamilies(VkSurfaceKHR surface);
        void createLogicalDevice();
        void createVmaAllocator();
    };

} // namespace Lizeral