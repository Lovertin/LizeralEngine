// source/runtime/function/render/rhi/vulkan/VulkanDevice.cpp

// 【魔法宏】：整个工程中只能有一个 cpp 文件定义它，它负责把 VMA 的实现展开
#define VMA_IMPLEMENTATION 
#include "VulkanDevice.h"
#include "VulkanContext.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>

namespace Lizeral {

    VulkanDevice::VulkanDevice(VulkanContext* context, VkSurfaceKHR surface) : m_context(context),m_surface(surface) {
        std::cout << "[VulkanDevice] Initializing Logical Device..." << std::endl;
        
        findQueueFamilies(surface);
        createLogicalDevice();
        createVmaAllocator();
    }

    VulkanDevice::~VulkanDevice() {
        std::cout << "[VulkanDevice] Destroying Logical Device..." << std::endl;
        
        if (m_allocator != VK_NULL_HANDLE) {
            vmaDestroyAllocator(m_allocator);
            m_allocator = VK_NULL_HANDLE;
        }

        if (m_logicalDevice != VK_NULL_HANDLE) {
            vkDestroyDevice(m_logicalDevice, nullptr);
            m_logicalDevice = VK_NULL_HANDLE;
        }
    }

    void VulkanDevice::WaitIdle() {
        if (m_logicalDevice != VK_NULL_HANDLE) {
            vkDeviceWaitIdle(m_logicalDevice);
        }
    }

    void VulkanDevice::findQueueFamilies(VkSurfaceKHR surface) {
        VkPhysicalDevice device = m_context->GetPhysicalDevice();
        uint32_t queueFamilyCount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, nullptr);
        std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
        vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, queueFamilies.data());

        for (uint32_t i = 0; i < queueFamilyCount; i++) {
            if (queueFamilies[i].queueFlags & VK_QUEUE_GRAPHICS_BIT) {
                m_queueIndices.graphicsFamily = i;
            }

            // ★ 检查该车间是否支持 Present (呈现)
            VkBool32 presentSupport = false;
            vkGetPhysicalDeviceSurfaceSupportKHR(device, i, surface, &presentSupport);
            if (presentSupport) {
                m_queueIndices.presentFamily = i;
            }

            if (m_queueIndices.isComplete()) break;
        }

        if (!m_queueIndices.isComplete()) throw std::runtime_error("Failed to find suitable queue families!");
        std::cout << "[VulkanDevice] Graphics Queue: " << m_queueIndices.graphicsFamily.value() 
                  << ", Present Queue: " << m_queueIndices.presentFamily.value() << std::endl;
    }

    void VulkanDevice::createLogicalDevice() {
        std::vector<VkDeviceQueueCreateInfo> queueCreateInfos;
        // ★ 新增：使用 set 去重。如果图形车间和呈现车间是同一个（通常都是 Index 0），我们只需要创建一个队列
        std::set<uint32_t> uniqueQueueFamilies = { m_queueIndices.graphicsFamily.value(), m_queueIndices.presentFamily.value() };

        float queuePriority = 1.0f; 
        for (uint32_t queueFamily : uniqueQueueFamilies) {
            VkDeviceQueueCreateInfo queueCreateInfo{};
            queueCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
            queueCreateInfo.queueFamilyIndex = queueFamily;
            queueCreateInfo.queueCount = 1;
            queueCreateInfo.pQueuePriorities = &queuePriority;
            queueCreateInfos.push_back(queueCreateInfo);
        }

        VkPhysicalDeviceFeatures deviceFeatures{};
        deviceFeatures.samplerAnisotropy = VK_TRUE; 

        VkPhysicalDeviceVulkan12Features features12{};
        features12.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES;
        features12.bufferDeviceAddress = VK_TRUE; 
        features12.scalarBlockLayout = VK_TRUE;

        VkPhysicalDeviceVulkan13Features features13{};
        features13.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES;
        features13.maintenance4 = VK_TRUE;

        // 新增：2. 开启 Vulkan 扩展的 Mesh Shader 特性
        VkPhysicalDeviceMeshShaderFeaturesEXT meshFeatures{};
        meshFeatures.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_MESH_SHADER_FEATURES_EXT;
        meshFeatures.meshShader = VK_TRUE; // 开启 Mesh Shader
        meshFeatures.taskShader = VK_TRUE; // 顺便把 Task Shader 也开了，以后剔除用
        
        // 核心魔法：把特性像铁链一样串起来传给 Vulkan (pNext 链)
        features13.pNext = &meshFeatures;
        features12.pNext = &meshFeatures;  // 把 meshFeatures 挂在 features12 后面

        VkDeviceCreateInfo createInfo{};
        createInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
        createInfo.queueCreateInfoCount = static_cast<uint32_t>(queueCreateInfos.size());
        createInfo.pQueueCreateInfos = queueCreateInfos.data();
        createInfo.pEnabledFeatures = &deviceFeatures;

        createInfo.pNext = &features12; 


        // ★ 新增：为了在窗口上画图，必须开启 Swapchain (交换链) 扩展！
        const std::vector<const char*> deviceExtensions = { 
            VK_KHR_SWAPCHAIN_EXTENSION_NAME ,
            VK_EXT_MESH_SHADER_EXTENSION_NAME
        };
        createInfo.enabledExtensionCount = static_cast<uint32_t>(deviceExtensions.size());
        createInfo.ppEnabledExtensionNames = deviceExtensions.data();

        if (vkCreateDevice(m_context->GetPhysicalDevice(), &createInfo, nullptr, &m_logicalDevice) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create logical Vulkan device!");
        }

        // 获取队列把柄
        vkGetDeviceQueue(m_logicalDevice, m_queueIndices.graphicsFamily.value(), 0, &m_graphicsQueue);
        // ★ 新增：获取呈现队列把柄
        vkGetDeviceQueue(m_logicalDevice, m_queueIndices.presentFamily.value(), 0, &m_presentQueue);

        std::cout << "[VulkanDevice] Logical Device created successfully!" << std::endl;
    }

    void VulkanDevice::createVmaAllocator() {
        VmaAllocatorCreateInfo allocatorInfo = {};
        allocatorInfo.physicalDevice = m_context->GetPhysicalDevice();
        allocatorInfo.device = m_logicalDevice;
        allocatorInfo.instance = (VkInstance)m_context->GetNativeInstance();
        allocatorInfo.vulkanApiVersion = VK_API_VERSION_1_3; // 保持与 Context 统一
        allocatorInfo.flags = VMA_ALLOCATOR_CREATE_BUFFER_DEVICE_ADDRESS_BIT;

        if (vmaCreateAllocator(&allocatorInfo, &m_allocator) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create VMA Allocator!");
        }

        std::cout << "[VulkanDevice] VMA Allocator initialized successfully!" << std::endl;
    }

} // namespace Lizeral