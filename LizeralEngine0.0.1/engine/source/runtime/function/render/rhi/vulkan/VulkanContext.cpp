#include "VulkanContext.h"
#include <iostream>
#include <stdexcept>
#include <cstring>

namespace Lizeral {

    // ========================================================
    // 宏配置：Debug 模式开启验证层，Release 模式关闭以追求极致性能
    // ========================================================
    const std::vector<const char*> validationLayers = {
        "VK_LAYER_KHRONOS_validation"
    };

#ifdef NDEBUG
    const bool enableValidationLayers = false;
#else
    const bool enableValidationLayers = true;
#endif

    // ========================================================
    // Vulkan 官方的 Debug 回调代理函数（因为它是扩展函数，需手动加载函数指针）
    // ========================================================
    VkResult CreateDebugUtilsMessengerEXT(VkInstance instance, const VkDebugUtilsMessengerCreateInfoEXT* pCreateInfo, const VkAllocationCallbacks* pAllocator, VkDebugUtilsMessengerEXT* pDebugMessenger) {
        auto func = (PFN_vkCreateDebugUtilsMessengerEXT) vkGetInstanceProcAddr(instance, "vkCreateDebugUtilsMessengerEXT");
        if (func != nullptr) {
            return func(instance, pCreateInfo, pAllocator, pDebugMessenger);
        } else {
            return VK_ERROR_EXTENSION_NOT_PRESENT;
        }
    }

    void DestroyDebugUtilsMessengerEXT(VkInstance instance, VkDebugUtilsMessengerEXT debugMessenger, const VkAllocationCallbacks* pAllocator) {
        auto func = (PFN_vkDestroyDebugUtilsMessengerEXT) vkGetInstanceProcAddr(instance, "vkDestroyDebugUtilsMessengerEXT");
        if (func != nullptr) {
            func(instance, debugMessenger, pAllocator);
        }
    }

    // 当 Vulkan 发生错误时，会触发这个回调函数打印红字！
    static VKAPI_ATTR VkBool32 VKAPI_CALL debugCallback(
        VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity,
        VkDebugUtilsMessageTypeFlagsEXT messageType,
        const VkDebugUtilsMessengerCallbackDataEXT* pCallbackData,
        void* pUserData) {
        
        std::cerr << "\n[Vulkan Validation Layer] " << pCallbackData->pMessage << "\n" << std::endl;
        return VK_FALSE;
    }

    // ========================================================
    // Context 实现部分
    // ========================================================

    void VulkanContext::Initialize(const std::string& appName, const std::vector<const char*>& windowExtensions) {
        std::cout << "[VulkanContext] Initializing Vulkan API..." << std::endl;
        
        // ★ 直接把 windowExtensions 传给 createInstance
        createInstance(appName, windowExtensions); 
        
        setupDebugMessenger();
        pickPhysicalDevice(); // 自动寻找 RTX 4060！
    }

    void VulkanContext::Shutdown() {
        if (enableValidationLayers) {
            DestroyDebugUtilsMessengerEXT(m_instance, m_debugMessenger, nullptr);
        }
        if (m_instance != VK_NULL_HANDLE) {
            std::cout << "[VulkanContext] Destroying Vulkan Instance..." << std::endl;
            vkDestroyInstance(m_instance, nullptr);
            m_instance = VK_NULL_HANDLE;
        }
    }

    void VulkanContext::createInstance(const std::string& appName, const std::vector<const char*>& windowExtensions) {
        if (enableValidationLayers && !checkValidationLayerSupport()) {
            throw std::runtime_error("Validation layers requested, but not available!");
        }

        VkApplicationInfo appInfo{};
        appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
        appInfo.pApplicationName = appName.c_str();
        appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
        appInfo.pEngineName = "Lizeral Engine";
        appInfo.engineVersion = VK_MAKE_VERSION(1, 0, 0);
        appInfo.apiVersion = VK_API_VERSION_1_3; 

        VkInstanceCreateInfo createInfo{};
        createInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
        createInfo.pApplicationInfo = &appInfo;

        // ★ 核心修复：在这里获取基础扩展，并合并窗口扩展！
        auto extensions = getRequiredExtensions();
        extensions.insert(extensions.end(), windowExtensions.begin(), windowExtensions.end());

        // 打印出来验明正身！
        std::cout << "[VulkanContext] Enabling " << extensions.size() << " Instance Extensions:" << std::endl;
        for (const auto& ext : extensions) {
            std::cout << "  - " << ext << std::endl;
        }

        createInfo.enabledExtensionCount = static_cast<uint32_t>(extensions.size());
        createInfo.ppEnabledExtensionNames = extensions.data();

        // 挂载验证层
        if (enableValidationLayers) {
            createInfo.enabledLayerCount = static_cast<uint32_t>(validationLayers.size());
            createInfo.ppEnabledLayerNames = validationLayers.data();
            std::cout << "[VulkanContext] Validation Layers ENABLED." << std::endl;
        } else {
            createInfo.enabledLayerCount = 0;
        }

        if (vkCreateInstance(&createInfo, nullptr, &m_instance) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create Vulkan instance!");
        }
        std::cout << "[VulkanContext] Successfully created Vulkan Instance!" << std::endl;
    }

    void VulkanContext::setupDebugMessenger() {
        if (!enableValidationLayers) return;

        VkDebugUtilsMessengerCreateInfoEXT createInfo{};
        createInfo.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT;
        // 抓取 Warning 和 Error 级别的报错
        createInfo.messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
        createInfo.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
        createInfo.pfnUserCallback = debugCallback;

        if (CreateDebugUtilsMessengerEXT(m_instance, &createInfo, nullptr, &m_debugMessenger) != VK_SUCCESS) {
            throw std::runtime_error("Failed to set up debug messenger!");
        }
        std::cout << "[VulkanContext] Debug Messenger setup successful." << std::endl;
    }

    void VulkanContext::pickPhysicalDevice() {
        uint32_t deviceCount = 0;
        vkEnumeratePhysicalDevices(m_instance, &deviceCount, nullptr);
        if (deviceCount == 0) throw std::runtime_error("Failed to find GPUs with Vulkan support!");

        std::vector<VkPhysicalDevice> devices(deviceCount);
        vkEnumeratePhysicalDevices(m_instance, &deviceCount, devices.data());

        std::cout << "[VulkanContext] Scanning " << deviceCount << " Physical Devices:" << std::endl;
        
        // 遍历所有显卡，优先选择独立显卡
        for (const auto& device : devices) {
            VkPhysicalDeviceProperties properties;
            vkGetPhysicalDeviceProperties(device, &properties);
            
            std::cout << "  - Found GPU: " << properties.deviceName << " (Type: " << properties.deviceType << ")" << std::endl;

            if (properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) {
                m_physicalDevice = device;
                std::cout << "[VulkanContext] >>> Selected Discrete GPU: " << properties.deviceName << " <<<" << std::endl;
                break; // 找到独显，直接选定并跳出
            }
        }

        // 如果没有独显，回退到第一个设备（通常是核显）
        if (m_physicalDevice == VK_NULL_HANDLE) {
            m_physicalDevice = devices[0];
            VkPhysicalDeviceProperties properties;
            vkGetPhysicalDeviceProperties(m_physicalDevice, &properties);
            std::cout << "[VulkanContext] >>> Warning: Fallback to GPU: " << properties.deviceName << " <<<" << std::endl;
        }
    }

    bool VulkanContext::checkValidationLayerSupport() {
        uint32_t layerCount;
        vkEnumerateInstanceLayerProperties(&layerCount, nullptr);
        std::vector<VkLayerProperties> availableLayers(layerCount);
        vkEnumerateInstanceLayerProperties(&layerCount, availableLayers.data());

        for (const char* layerName : validationLayers) {
            bool layerFound = false;
            for (const auto& layerProperties : availableLayers) {
                if (strcmp(layerName, layerProperties.layerName) == 0) {
                    layerFound = true;
                    break;
                }
            }
            if (!layerFound) return false;
        }
        return true;
    }

    std::vector<const char*> VulkanContext::getRequiredExtensions() {
        std::vector<const char*> extensions;
        if (enableValidationLayers) {
            extensions.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
        }
        return extensions;
    }

} // namespace Lizeral