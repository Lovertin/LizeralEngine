#pragma once
#include "runtime/function/render/rhi/RHIContext.h"
#include <vulkan/vulkan.h>
#include <vector>
#include <string>

namespace Lizeral {

    class VulkanContext : public IRHIContext {
    public:
        VulkanContext() = default;
        ~VulkanContext() override { Shutdown(); }

        void Initialize(const std::string& appName, const std::vector<const char*>& windowExtensions = {});
        void Shutdown() override;
        
        void* GetNativeInstance() const override { return m_instance; }

        VkPhysicalDevice GetPhysicalDevice() const { return m_physicalDevice; }

    private:
        VkInstance m_instance { VK_NULL_HANDLE };
        VkDebugUtilsMessengerEXT m_debugMessenger { VK_NULL_HANDLE };
        VkPhysicalDevice m_physicalDevice { VK_NULL_HANDLE }; 

        void createInstance(const std::string& appName,const std::vector<const char*>& windowExtensions);
        void setupDebugMessenger();
        void pickPhysicalDevice();

        bool checkValidationLayerSupport();
        std::vector<const char*> getRequiredExtensions();
        void printAvailableExtensions();
    };

} // namespace Lizeral