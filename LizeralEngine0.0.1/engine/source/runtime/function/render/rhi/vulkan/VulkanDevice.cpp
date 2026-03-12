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
        // ★ 使用 set 去重。如果图形车间和呈现车间是同一个（通常都是 Index 0），我们只需要创建一个队列
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
        deviceFeatures.shaderInt64 = VK_TRUE;

        // ========================================================
        // ★ 0. 新增：光线追踪相关特性 (Ray Query & AS)
        // ========================================================
        VkPhysicalDeviceRayQueryFeaturesKHR rayQueryFeatures{};
        rayQueryFeatures.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_RAY_QUERY_FEATURES_KHR;
        rayQueryFeatures.rayQuery = VK_TRUE; // 开启在 Shader 中发射射线的能力

        VkPhysicalDeviceAccelerationStructureFeaturesKHR asFeatures{};
        asFeatures.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_FEATURES_KHR;
        asFeatures.accelerationStructure = VK_TRUE; // 开启构建加速结构的能力
        asFeatures.pNext = &rayQueryFeatures;       // 链条：AS -> RayQuery

        // 1. Mesh Shader 特性
        VkPhysicalDeviceMeshShaderFeaturesEXT meshFeatures{};
        meshFeatures.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_MESH_SHADER_FEATURES_EXT;
        meshFeatures.meshShader = VK_TRUE; 
        meshFeatures.taskShader = VK_TRUE; 
        meshFeatures.pNext = &asFeatures;           // ★ 链条：Mesh -> AS

        // 2. Vulkan 1.3 特性
        VkPhysicalDeviceVulkan13Features features13{};
        features13.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES;
        features13.dynamicRendering = VK_TRUE; 
        features13.synchronization2 = VK_TRUE; 
        features13.maintenance4 = VK_TRUE;
        features13.pNext = &meshFeatures;           // 链条第 1 环 -> Mesh

        // 3. Vulkan 1.2 特性
        VkPhysicalDeviceVulkan12Features features12{};
        features12.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES;
        features12.bufferDeviceAddress = VK_TRUE;  // ★ 光追极度依赖 BDA，你已经开启了，很棒！
        features12.scalarBlockLayout = VK_TRUE;
        features12.descriptorIndexing = VK_TRUE;                              
        features12.shaderSampledImageArrayNonUniformIndexing = VK_TRUE;       
        features12.descriptorBindingPartiallyBound = VK_TRUE;                 
        features12.runtimeDescriptorArray = VK_TRUE;
        features12.pNext = &features13;            // 链条第 2 环 -> features13

        // 4. Vulkan 1.1 特性 (开启 Multiview)
        VkPhysicalDeviceVulkan11Features features11{};
        features11.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES;
        features11.multiview = VK_TRUE;
        features11.pNext = &features12;            // 链条第 3 环 -> features12

        VkDeviceCreateInfo createInfo{};
        createInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
        createInfo.queueCreateInfoCount = static_cast<uint32_t>(queueCreateInfos.size());
        createInfo.pQueueCreateInfos = queueCreateInfos.data();
        createInfo.pEnabledFeatures = &deviceFeatures;
        
        // 最终提交链条的头部改为 features11！
        createInfo.pNext = &features11;

        // ★ 注意：这里帮你补上了 MESH_SHADER 扩展后面的逗号！
        const std::vector<const char*> deviceExtensions = { 
            VK_KHR_SWAPCHAIN_EXTENSION_NAME,
            VK_EXT_MESH_SHADER_EXTENSION_NAME,             // <-- 刚才漏了这里的逗号
            VK_KHR_ACCELERATION_STRUCTURE_EXTENSION_NAME,  // 加速结构
            VK_KHR_RAY_QUERY_EXTENSION_NAME,               // Ray Query
            VK_KHR_DEFERRED_HOST_OPERATIONS_EXTENSION_NAME // AS 构建依赖它
        };
        createInfo.enabledExtensionCount = static_cast<uint32_t>(deviceExtensions.size());
        createInfo.ppEnabledExtensionNames = deviceExtensions.data();

        if (vkCreateDevice(m_context->GetPhysicalDevice(), &createInfo, nullptr, &m_logicalDevice) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create logical Vulkan device!");
        }

        // 获取队列把柄
        vkGetDeviceQueue(m_logicalDevice, m_queueIndices.graphicsFamily.value(), 0, &m_graphicsQueue);
        // 获取呈现队列把柄
        vkGetDeviceQueue(m_logicalDevice, m_queueIndices.presentFamily.value(), 0, &m_presentQueue);

        std::cout << "[VulkanDevice] Logical Device created successfully with Ray Tracing enabled!" << std::endl;
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