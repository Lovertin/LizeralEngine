#include <iostream>
#include <exception>
#include <fstream>
#include <vector>
#include <cstring>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"
#include "runtime/function/render/rhi/vulkan/VulkanPipelineBuilder.h"

// 引入我们的 Meshlet 兵工厂
#include "runtime/function/render/MeshletBuilder/MeshletModelBuilder.h" 

const uint32_t WIDTH = 1280;
const uint32_t HEIGHT = 720;

// ★ 请根据你的实际路径修改这里！
#define SHADER_DIR "C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/shader/"
#define MODEL_PATH "C:/Lizeral Engine/LizeralEngine0.0.1/asset/mazda_glb.glb"

// ====================================================================
// 泛型 BDA 分配器
// ====================================================================
template <typename T>
uint64_t createBDABuffer(VkDevice device, VkPhysicalDevice physicalDevice, const std::vector<T>& data, VkBuffer& outBuffer, VkDeviceMemory& outMemory) {
    VkBufferCreateInfo bufferInfo{};
    bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bufferInfo.size = data.size() * sizeof(T);
    bufferInfo.usage = VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
    bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    vkCreateBuffer(device, &bufferInfo, nullptr, &outBuffer);

    VkMemoryRequirements memRequirements;
    vkGetBufferMemoryRequirements(device, outBuffer, &memRequirements);

    VkPhysicalDeviceMemoryProperties memProperties;
    vkGetPhysicalDeviceMemoryProperties(physicalDevice, &memProperties);
    uint32_t memoryTypeIndex = 0;
    for (uint32_t i = 0; i < memProperties.memoryTypeCount; i++) {
        if ((memRequirements.memoryTypeBits & (1 << i)) && 
            (memProperties.memoryTypes[i].propertyFlags & (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)) == (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)) {
            memoryTypeIndex = i;
            break;
        }
    }

    VkMemoryAllocateFlagsInfo flagsInfo{};
    flagsInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO;
    flagsInfo.flags = VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT; 

    VkMemoryAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocInfo.allocationSize = memRequirements.size;
    allocInfo.memoryTypeIndex = memoryTypeIndex;
    allocInfo.pNext = &flagsInfo;

    vkAllocateMemory(device, &allocInfo, nullptr, &outMemory);
    vkBindBufferMemory(device, outBuffer, outMemory, 0);

    void* mappedData;
    vkMapMemory(device, outMemory, 0, bufferInfo.size, 0, &mappedData);
    memcpy(mappedData, data.data(), (size_t)bufferInfo.size);
    vkUnmapMemory(device, outMemory);

    VkBufferDeviceAddressInfo addressInfo{};
    addressInfo.sType = VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO;
    addressInfo.buffer = outBuffer;
    return vkGetBufferDeviceAddress(device, &addressInfo);
}

// ====================================================================
// GPU 任务信使
// ====================================================================
struct PushConstants {
    float transform[4];     // x, y, z 偏移量， w 为统一缩放
    uint64_t vertexBuffer;  // 全局顶点数组物理地址
    uint64_t meshletBuffer; // Meshlet 描述符数组物理地址
    uint64_t indexBuffer;   // 微型索引数组物理地址
};

std::vector<char> readFile(const std::string& filename) {
    std::ifstream file(filename, std::ios::ate | std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Failed to open file: " + filename);
    size_t fileSize = (size_t)file.tellg();
    std::vector<char> buffer(fileSize);
    file.seekg(0);
    file.read(buffer.data(), fileSize);
    file.close();
    return buffer;
}

VkShaderModule createShaderModule(VkDevice device, const std::vector<char>& code) {
    VkShaderModuleCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    createInfo.codeSize = code.size();
    createInfo.pCode = reinterpret_cast<const uint32_t*>(code.data());
    VkShaderModule shaderModule;
    if (vkCreateShaderModule(device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create shader module!");
    }
    return shaderModule;
}

int main() {
    GLFWwindow* window = nullptr;
    VkSurfaceKHR surface = VK_NULL_HANDLE;
    VkInstance instance = VK_NULL_HANDLE;

    try {
        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API); 
        glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
        window = glfwCreateWindow(WIDTH, HEIGHT, "Lizeral Engine - Dynamic Rendering & Meshlet", nullptr, nullptr);

        Lizeral::VulkanContext vulkanContext;
        uint32_t glfwExtensionCount = 0;
        const char** glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);
        std::vector<const char*> requiredExtensions(glfwExtensions, glfwExtensions + glfwExtensionCount);
        vulkanContext.Initialize("Lizeral Sandbox", requiredExtensions);
        instance = (VkInstance)vulkanContext.GetNativeInstance();

        if (glfwCreateWindowSurface(instance, window, nullptr, &surface) != VK_SUCCESS) throw std::runtime_error("Failed to create window surface!");

        {
            Lizeral::VulkanDevice vulkanDevice(&vulkanContext, surface);
            
            auto CmdDrawMeshTasksEXT = (PFN_vkCmdDrawMeshTasksEXT)vkGetDeviceProcAddr(vulkanDevice.GetNativeDevice(), "vkCmdDrawMeshTasksEXT");
            if (!CmdDrawMeshTasksEXT) {
                throw std::runtime_error("Failed to load vkCmdDrawMeshTasksEXT! Please ensure Vulkan 1.2+ and VK_EXT_mesh_shader are enabled.");
            }

            // ====================================================================
            // 1. 资产管线发力：加载模型并用 MeshOptimizer 切分
            // ====================================================================
            Lizeral::MeshletModelBuilder builder;
            if (!builder.LoadAndSliceModel(MODEL_PATH)) {
                throw std::runtime_error("Failed to load or slice model!");
            }

            // ====================================================================
            // 2. 分配 BDA 显存并装填通讯包
            // ====================================================================
            VkBuffer vBuf, mBuf, iBuf;
            VkDeviceMemory vMem, mMem, iMem;
            
            uint64_t vAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetVertices(), vBuf, vMem);
            uint64_t mAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetMeshlets(), mBuf, mMem);
            uint64_t iAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetMicroIndices(), iBuf, iMem);

            PushConstants pushData{};
            pushData.transform[0] = 0.0f;    
            pushData.transform[1] = 0.0f;    
            pushData.transform[2] = 0.0f;    
            pushData.transform[3] = 10.5f;  // ★ 之前测试出能完美看见侧脸的缩放系数
            pushData.vertexBuffer = vAddr;
            pushData.meshletBuffer = mAddr;
            pushData.indexBuffer = iAddr;
            
            Lizeral::VulkanRenderer renderer(&vulkanContext, &vulkanDevice, window);

            auto meshShaderCode = readFile(std::string(SHADER_DIR) + "car_mesh.spv"); 
            auto fragShaderCode = readFile(std::string(SHADER_DIR) + "car_frag.spv");
            VkShaderModule meshShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), meshShaderCode);
            VkShaderModule fragShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), fragShaderCode);

            VkPushConstantRange pushConstantRange{};
            pushConstantRange.stageFlags = VK_SHADER_STAGE_MESH_BIT_EXT;
            pushConstantRange.offset = 0;
            pushConstantRange.size = sizeof(PushConstants); 

            VkPipelineLayoutCreateInfo pipelineLayoutInfo{};
            pipelineLayoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
            pipelineLayoutInfo.pushConstantRangeCount = 1;
            pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;
            VkPipelineLayout pipelineLayout;
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &pipelineLayoutInfo, nullptr, &pipelineLayout);

            Lizeral::VulkanPipelineBuilder builderPipeline;
            
            // ★ 解放的管线：不传 RenderPass，只传两个格式！
            VkPipeline graphicsPipeline = builderPipeline
                .AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShaderModule) 
                .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragShaderModule)
                .SetPipelineLayout(pipelineLayout)
                .Build(&vulkanDevice, renderer.GetSwapchainFormat(), VK_FORMAT_D32_SFLOAT);

            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), fragShaderModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), meshShaderModule, nullptr);

            std::cout << "\n[Sandbox] Entering Dynamic Rendering Loop! SIMT Power UNLEASHED!" << std::endl;
            
            // ====================================================================
            // 4. 次世代动态渲染主循环
            // ====================================================================
            while (!glfwWindowShouldClose(window)) {
                glfwPollEvents();

                VkCommandBuffer cmd = renderer.BeginFrame();
                
                if (cmd != VK_NULL_HANDLE) {
                    
                    // ★ 动态渲染开启！随叫随到，指哪画哪！
                    renderer.BeginRendering(cmd);
                    
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, graphicsPipeline);
                    vkCmdPushConstants(cmd, pipelineLayout, VK_SHADER_STAGE_MESH_BIT_EXT, 0, sizeof(PushConstants), &pushData);
                    
                    // 极致并发发射
                    CmdDrawMeshTasksEXT(cmd, builder.GetMeshlets().size(), 1, 1);

                    // ★ 动态渲染结束！
                    renderer.EndRendering(cmd);

                    renderer.EndFrame();
                }
            }
            
            vkDeviceWaitIdle(vulkanDevice.GetNativeDevice());
            
            // ====================================================================
            // 5. 资源清理
            // ====================================================================
            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), graphicsPipeline, nullptr);
            vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), pipelineLayout, nullptr);
            
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), vBuf, nullptr);
            vkFreeMemory(vulkanDevice.GetNativeDevice(), vMem, nullptr);
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), mBuf, nullptr);
            vkFreeMemory(vulkanDevice.GetNativeDevice(), mMem, nullptr);
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), iBuf, nullptr);
            vkFreeMemory(vulkanDevice.GetNativeDevice(), iMem, nullptr);
        } 

        vkDestroySurfaceKHR(instance, surface, nullptr);
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL ERROR] " << e.what() << std::endl;
    }

    if (window) {
        glfwDestroyWindow(window);
        glfwTerminate();
    }
    return 0;
}