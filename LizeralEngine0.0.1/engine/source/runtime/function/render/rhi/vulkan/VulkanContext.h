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

        // 获取选中的显卡（后续创建 Device 时要用到）
        VkPhysicalDevice GetPhysicalDevice() const { return m_physicalDevice; }

    private:
        VkInstance m_instance { VK_NULL_HANDLE };
        VkDebugUtilsMessengerEXT m_debugMessenger { VK_NULL_HANDLE };
        VkPhysicalDevice m_physicalDevice { VK_NULL_HANDLE }; // 我们选中的显卡

        // 内部流程函数
        void createInstance(const std::string& appName,const std::vector<const char*>& windowExtensions);
        void setupDebugMessenger();
        void pickPhysicalDevice();

        // 辅助检查函数
        bool checkValidationLayerSupport();
        std::vector<const char*> getRequiredExtensions();
        void printAvailableExtensions();
    };

} // namespace Lizeral