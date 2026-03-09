#include <iostream>
#include <exception>
#include <fstream>
#include <vector>
#include <cstring>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h" // ★ 引入我们的渲染调度器
#include "runtime/function/render/rhi/vulkan/VulkanPipelineBuilder.h"
#include "runtime/function/render/rhi/vulkan/VulkanRenderPass.h"

const uint32_t WIDTH = 1280;
const uint32_t HEIGHT = 720;

// ★ 请把这里换成你编译出的 .spv 文件所在的文件夹绝对路径
#define SHADER_DIR "C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/shader/"

// ====================================================================
// 顶点结构体与 BDA 显存分配器 (保持不变)
// ====================================================================
struct Vertex {
    float pos[3];
    float color[3];
};

uint64_t createBDAVertexBuffer(VkDevice device, VkPhysicalDevice physicalDevice, const std::vector<Vertex>& vertices, VkBuffer& outBuffer, VkDeviceMemory& outMemory) {
    VkBufferCreateInfo bufferInfo{};
    bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bufferInfo.size = vertices.size() * sizeof(Vertex);
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

    void* data;
    vkMapMemory(device, outMemory, 0, bufferInfo.size, 0, &data);
    memcpy(data, vertices.data(), (size_t)bufferInfo.size);
    vkUnmapMemory(device, outMemory);

    VkBufferDeviceAddressInfo addressInfo{};
    addressInfo.sType = VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO;
    addressInfo.buffer = outBuffer;
    return vkGetBufferDeviceAddress(device, &addressInfo);
}

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
        glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE); // ★ 设为 TRUE，测试我们渲染器的自动缩放功能！
        window = glfwCreateWindow(WIDTH, HEIGHT, "Lizeral Engine - Renderer Powered", nullptr, nullptr);

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
                throw std::runtime_error("Failed to load vkCmdDrawMeshTasksEXT!");
            }

            // ====================================================================
            // 1. 准备 BDA 顶点数据
            // ====================================================================
            std::vector<Vertex> triangleVertices = {
                {{ 0.0f, -0.5f, 0.0f}, {1.0f, 0.0f, 0.0f}},
                {{ 0.5f,  0.5f, 0.0f}, {0.0f, 1.0f, 0.0f}},
                {{-0.5f,  0.5f, 0.0f}, {0.0f, 0.0f, 1.0f}} 
            };
            VkBuffer vertexBuffer;
            VkDeviceMemory vertexMemory;
            uint64_t bdaPointer = createBDAVertexBuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), triangleVertices, vertexBuffer, vertexMemory);

            // ====================================================================
            // 2. 引擎系统启动：创建渲染器！
            // ====================================================================
            Lizeral::VulkanRenderer renderer(&vulkanContext, &vulkanDevice, window);

            // ====================================================================
            // 3. 准备管线
            // ====================================================================
            auto meshShaderCode = readFile(std::string(SHADER_DIR) + "triangle_mesh.spv"); 
            auto fragShaderCode = readFile(std::string(SHADER_DIR) + "triangle_frag.spv");
            VkShaderModule meshShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), meshShaderCode);
            VkShaderModule fragShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), fragShaderCode);

            VkPushConstantRange pushConstantRange{};
            pushConstantRange.stageFlags = VK_SHADER_STAGE_MESH_BIT_EXT;
            pushConstantRange.offset = 0;
            pushConstantRange.size = sizeof(uint64_t); 

            VkPipelineLayoutCreateInfo pipelineLayoutInfo{};
            pipelineLayoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
            pipelineLayoutInfo.pushConstantRangeCount = 1;
            pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;
            VkPipelineLayout pipelineLayout;
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &pipelineLayoutInfo, nullptr, &pipelineLayout);

            Lizeral::VulkanPipelineBuilder builder;
            VkPipeline graphicsPipeline = builder
                .AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShaderModule) 
                .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragShaderModule)
                .SetPipelineLayout(pipelineLayout)
                .Build(&vulkanDevice, renderer.GetSwapChainRenderPass()->GetNativeRenderPass()); // ★ 直接向 Renderer 要 RenderPass！

            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), fragShaderModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), meshShaderModule, nullptr);

            std::cout << "\n[Sandbox] Entering Render Loop! Mesh Shader Activated!" << std::endl;
            
            // ====================================================================
            // ★ 4. 极致清爽的现代引擎主循环
            // ====================================================================
            while (!glfwWindowShouldClose(window)) {
                glfwPollEvents();

                // 【步骤 A】：请求开启一帧 (底层帮你搞定同步锁、获取图片等所有脏活累活)
                VkCommandBuffer cmd = renderer.BeginFrame();
                
                // 如果窗口被最小化，BeginFrame 会返回 NULL，此时跳过绘制
                if (cmd != VK_NULL_HANDLE) {
                    
                    // 【步骤 B】：开启最终上屏的 RenderPass
                    renderer.BeginSwapChainRenderPass(cmd);
                    
                    // ============================================
                    // ★ 这里就是未来填充游戏渲染逻辑 (ECS遍历) 的地方！
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, graphicsPipeline);
                    vkCmdPushConstants(cmd, pipelineLayout, VK_SHADER_STAGE_MESH_BIT_EXT, 0, sizeof(uint64_t), &bdaPointer);
                    CmdDrawMeshTasksEXT(cmd, 1, 1, 1);
                    // ============================================

                    // 【步骤 C】：结束 RenderPass
                    renderer.EndSwapChainRenderPass(cmd);

                    // 【步骤 D】：请求结束一帧 (底层帮你提交 Queue 和 Present)
                    renderer.EndFrame();
                }
            }
            
            // ====================================================================
            // 5. 优雅退出与清理
            // ====================================================================
            vkDeviceWaitIdle(vulkanDevice.GetNativeDevice());
            
            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), graphicsPipeline, nullptr);
            vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), pipelineLayout, nullptr);
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), vertexBuffer, nullptr);
            vkFreeMemory(vulkanDevice.GetNativeDevice(), vertexMemory, nullptr);
            
        } // Renderer 等底层核心对象在此被自动安全释放

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