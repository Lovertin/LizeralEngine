#include <iostream>
#include <exception>
#include <fstream>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <memory>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

// --- Vulkan & Engine RHI ---
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"
#include "runtime/function/render/rhi/vulkan/VulkanPipelineBuilder.h"
#include "runtime/function/render/MeshletBuilder/MeshletModelBuilder.h" 
#include "runtime/function/render/rhi/vulkan/VulkanCommandPool.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandBuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanTexture.h"
#include "runtime/function/render/rhi/vulkan/VulkanDescriptorBuilder.h"
#include "runtime/function/Vulkan_res_type/VulkanModelResource.h"
#include "runtime/function/render/rhi/vulkan/VulkanBLAS.h"
#include "runtime/function/render/rhi/vulkan/VulkanTLAS.h"

// --- ECS & Camera System ---
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/ecs/components/Model/VulkanModelComponent.h" // ★ 引入新的哑巴数据组件
#include "runtime/function/input/input.h"
#include "runtime/function/render/CameraControlSystem/CameraControlSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h" 
#include "runtime/core/math/matrix4.h"

using namespace Lizeral;

// ====================================================================
// 1. 核心定义区
// ====================================================================
const uint32_t WIDTH = 1280;
const uint32_t HEIGHT = 720;
const std::string SHADER_DIR = "C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/shader/";

struct PushConstants {
    Matrix4x4 mvp; 
    Matrix4x4 model;
    Matrix4x4 prevMvp;
    uint64_t vertexBuffer;
    uint64_t meshletBuffer;
    uint64_t indexBuffer;
    uint64_t boundsBuffer; 
    uint64_t materialBuffer;
    uint32_t totalMeshlets;
    uint32_t textureOffset;
    Vector2 jitter;
};

struct LightingPushConstants {
    Lizeral::Matrix4x4 invViewProj; 
    Lizeral::Matrix4x4 viewProj;
    Lizeral::Vector3 cameraPos;     
    float padding;                  
};

struct GBufferAttachment {
    VkImage image = VK_NULL_HANDLE;
    VmaAllocation allocation = VK_NULL_HANDLE;
    VkImageView view = VK_NULL_HANDLE;
    VkFormat format;
};

// 用于统一清理 BDA 物理内存的辅助结构
struct BufferAllocation {
    VkBuffer buffer;
    VkDeviceMemory memory;
};

// ====================================================================
// 2. 工具函数与回调区
// ====================================================================
void glfwKeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    Key lizKey; bool found = true;
    switch (key) {
        case GLFW_KEY_W: lizKey = Key::W; break;
        case GLFW_KEY_A: lizKey = Key::A; break;
        case GLFW_KEY_S: lizKey = Key::S; break;
        case GLFW_KEY_D: lizKey = Key::D; break;
        case GLFW_KEY_Q: lizKey = Key::Q; break;
        case GLFW_KEY_E: lizKey = Key::E; break;
        case GLFW_KEY_LEFT_SHIFT: lizKey = Key::LEFT_SHIFT; break;
        case GLFW_KEY_ESCAPE: lizKey = Key::ESC; break;
        default: found = false; break;
    }
    if (found) Input::GetInstance().SetKeyDown(lizKey, action != GLFW_RELEASE);
}

void glfwCursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    Input::GetInstance().SetMousePosition(static_cast<float>(xpos), static_cast<float>(ypos));
}

void glfwMouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        Input::GetInstance().SetMouseButtonDown(MouseButton::Right, action != GLFW_RELEASE);
    }
}

std::vector<char> readFile(const std::string& filename) {
    std::ifstream file(filename, std::ios::ate | std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Failed to open file: " + filename);
    size_t fileSize = (size_t)file.tellg();
    std::vector<char> buffer(fileSize);
    file.seekg(0); file.read(buffer.data(), fileSize); file.close();
    return buffer;
}

VkShaderModule createShaderModule(VkDevice device, const std::vector<char>& code) {
    VkShaderModuleCreateInfo createInfo{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    createInfo.codeSize = code.size();
    createInfo.pCode = reinterpret_cast<const uint32_t*>(code.data());
    VkShaderModule shaderModule;
    if (vkCreateShaderModule(device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) 
        throw std::runtime_error("Failed to create shader module!");
    return shaderModule;
}

template <typename T>
uint64_t createBDABuffer(VkDevice device, VkPhysicalDevice physicalDevice, const std::vector<T>& data, VkBuffer& outBuffer, VkDeviceMemory& outMemory) {
    VkBufferCreateInfo bufferInfo{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bufferInfo.size = data.size() * sizeof(T);
    bufferInfo.usage = VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT | VK_BUFFER_USAGE_ACCELERATION_STRUCTURE_BUILD_INPUT_READ_ONLY_BIT_KHR;
    bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    vkCreateBuffer(device, &bufferInfo, nullptr, &outBuffer);

    VkMemoryRequirements memReqs;
    vkGetBufferMemoryRequirements(device, outBuffer, &memReqs);

    VkPhysicalDeviceMemoryProperties memProps;
    vkGetPhysicalDeviceMemoryProperties(physicalDevice, &memProps);
    uint32_t memoryTypeIndex = 0;
    for (uint32_t i = 0; i < memProps.memoryTypeCount; i++) {
        if ((memReqs.memoryTypeBits & (1 << i)) && 
            (memProps.memoryTypes[i].propertyFlags & (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)) == (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)) {
            memoryTypeIndex = i; break;
        }
    }

    VkMemoryAllocateFlagsInfo flagsInfo{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO};
    flagsInfo.flags = VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT; 
    VkMemoryAllocateInfo allocInfo{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    allocInfo.allocationSize = memReqs.size;
    allocInfo.memoryTypeIndex = memoryTypeIndex;
    allocInfo.pNext = &flagsInfo;

    vkAllocateMemory(device, &allocInfo, nullptr, &outMemory);
    vkBindBufferMemory(device, outBuffer, outMemory, 0);

    void* mappedData;
    vkMapMemory(device, outMemory, 0, bufferInfo.size, 0, &mappedData);
    memcpy(mappedData, data.data(), (size_t)bufferInfo.size);
    vkUnmapMemory(device, outMemory);

    VkBufferDeviceAddressInfo addressInfo{VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO};
    addressInfo.buffer = outBuffer;
    return vkGetBufferDeviceAddress(device, &addressInfo);
}

// 完善后的布局转换函数 (支持 G-Buffer 循环闭环)
void transitionImageLayout(VkCommandBuffer cmd, VkImage image, VkImageLayout oldLayout, VkImageLayout newLayout, VkImageAspectFlags aspectMask) {
    VkImageMemoryBarrier barrier{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
    barrier.oldLayout = oldLayout;
    barrier.newLayout = newLayout;
    barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.image = image;
    barrier.subresourceRange.aspectMask = aspectMask;
    barrier.subresourceRange.baseMipLevel = 0;
    barrier.subresourceRange.levelCount = 1;
    barrier.subresourceRange.baseArrayLayer = 0;
    barrier.subresourceRange.layerCount = 1;

    VkPipelineStageFlags sourceStage;
    VkPipelineStageFlags destinationStage;

    if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL) {
        barrier.srcAccessMask = 0; 
        barrier.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
        sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT; 
        destinationStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    } 
    // ★ 新增：支持从未定义直接转为只读 (用于 TAA 历史帧初始化)
    else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
        barrier.srcAccessMask = 0;
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
        destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    }
    else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL) {
        barrier.srcAccessMask = 0; 
        barrier.dstAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
        sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT; 
        destinationStage = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT;
    } 
    else if (oldLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
        barrier.srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT; 
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        sourceStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT; 
        destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    } 
    else if (oldLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
        barrier.srcAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT; 
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        sourceStage = VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT; 
        destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    } 
    else if (oldLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL) {
        barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT; 
        barrier.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
        sourceStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT; 
        destinationStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    } 
    else if (oldLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL) {
        barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT;
        barrier.dstAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
        sourceStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
        // 深度测试发生在这里
        destinationStage = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT; 
    }
    else {
        // 如果再报错，这里能打印出到底是哪两个 Layout 没对上
        std::cerr << "Unsupported transition: " << oldLayout << " -> " << newLayout << std::endl;
        throw std::invalid_argument("Unsupported layout transition!");
    }

    vkCmdPipelineBarrier(cmd, sourceStage, destinationStage, 0, 0, nullptr, 0, nullptr, 1, &barrier);
}

float CreateHaltonSequence(uint32_t index, uint32_t base) {
    float f = 1.0f;
    float result = 0.0f;
    uint32_t i = index;
    while (i > 0) {
        f = f / static_cast<float>(base);
        result = result + f * static_cast<float>(i % base);
        i = static_cast<uint32_t>(std::floor(static_cast<float>(i) / static_cast<float>(base)));
    }
    return result;
}

// ====================================================================
// 3. 主程序入口
// ====================================================================
int main() {
    GLFWwindow* window = nullptr;
    VkSurfaceKHR surface = VK_NULL_HANDLE;
    VkInstance instance = VK_NULL_HANDLE;

    try {
        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API); 
        glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
        window = glfwCreateWindow(WIDTH, HEIGHT, "Lizeral Engine - Mesh Shader Rendering", nullptr, nullptr);

        glfwSetKeyCallback(window, glfwKeyCallback);
        glfwSetCursorPosCallback(window, glfwCursorPosCallback);
        glfwSetMouseButtonCallback(window, glfwMouseButtonCallback);

        VulkanContext vulkanContext;
        uint32_t glfwExtensionCount = 0;
        const char** glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);
        std::vector<const char*> requiredExtensions(glfwExtensions, glfwExtensions + glfwExtensionCount);
        vulkanContext.Initialize("Lizeral Sandbox", requiredExtensions);
        instance = (VkInstance)vulkanContext.GetNativeInstance();

        if (glfwCreateWindowSurface(instance, window, nullptr, &surface) != VK_SUCCESS) 
            throw std::runtime_error("Failed to create window surface!");

        {
            VulkanDevice vulkanDevice(&vulkanContext, surface);
            VulkanRenderer renderer(&vulkanContext, &vulkanDevice, window);
            Lizeral::VulkanCommandPool resourceCommandPool(&vulkanDevice);

            // =======================================================
            // ★ 核心架构：全局资源与贴图宇宙池
            // =======================================================
            std::unordered_map<std::string, Lizeral::VulkanModelResource> g_ModelCache;
            std::vector<std::unique_ptr<Lizeral::VulkanTexture>> g_GlobalTextures; 
            std::vector<VkDescriptorImageInfo> g_GlobalImageInfos;
            std::vector<BufferAllocation> g_AllocatedBuffers; // 统一记录并回收所有 BDA 缓冲

            VulkanTLAS tlas(&vulkanDevice);

            //加载模型的万用lambda函数
            auto getOrLoadModel = [&](const std::string& path) -> Lizeral::VulkanModelResource {
                if (g_ModelCache.find(path) != g_ModelCache.end()) return g_ModelCache[path];

                std::cout << "[RenderSystem] Loading new model to GPU: " << path << std::endl;
                
                uint32_t currentTexOffset = static_cast<uint32_t>(g_GlobalTextures.size());
                
                MeshletModelBuilder builder;
                if (!builder.LoadAndSliceModel(path, currentTexOffset)) {
                    throw std::runtime_error("Failed to load GLB model: " + path);
                }
                
                for (const auto& texData : builder.GetAllTextures()) {
                    auto texture = std::make_unique<Lizeral::VulkanTexture>(&vulkanDevice, &resourceCommandPool, texData.data(), texData.size());
                    
                    VkDescriptorImageInfo info{};
                    info.imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
                    info.imageView = texture->GetImageView();
                    info.sampler = texture->GetSampler(); 
                    g_GlobalImageInfos.push_back(info);
                    g_GlobalTextures.push_back(std::move(texture));
                }

                Lizeral::VulkanModelResource res;
                VkBuffer vBuf, mBuf, iBuf, bBuf, matBuf;
                VkDeviceMemory vMem, mMem, iMem, bMem, matMem;
                
                res.vAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetVertices(), vBuf, vMem);
                res.mAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetMeshlets(), mBuf, mMem);
                res.iAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetMicroIndices(), iBuf, iMem);
                res.bAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetBounds(), bBuf, bMem);
                res.matAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetMaterials(), matBuf, matMem);
                
                // 登记到销毁清单
                g_AllocatedBuffers.push_back({vBuf, vMem}); g_AllocatedBuffers.push_back({mBuf, mMem});
                g_AllocatedBuffers.push_back({iBuf, iMem}); g_AllocatedBuffers.push_back({bBuf, bMem});
                g_AllocatedBuffers.push_back({matBuf, matMem});

                res.totalMeshlets = builder.GetMeshlets().size();
                res.textureOffset = currentTexOffset;
                res.textureCount = builder.GetAllTextures().size();

                // ========================================================
                // ★ 新增：在这里构建 BLAS (Bottom-Level Acceleration Structure)
                // ========================================================
                const auto& vertices = builder.GetVertices();
                const auto& indices = builder.GetMicroIndices();
                
                uint32_t vertexCount = static_cast<uint32_t>(vertices.size());
                uint32_t indexCount = static_cast<uint32_t>(indices.size());
                
                // 自动推断单个顶点的字节大小 (Stride)
                uint32_t vertexStride = vertices.empty() ? 0 : static_cast<uint32_t>(sizeof(vertices[0]));

                if (vertexCount > 0 && indexCount > 0) {
                    std::cout << "[VulkanBLAS] Triggering BLAS build for: " << path << std::endl;
                    // 实例化并构建 BLAS，GPU 会在这里被命令记录并执行构建
                    res.blas = std::make_shared<VulkanBLAS>(
                        &vulkanDevice, 
                        &resourceCommandPool, 
                        res.vAddr, vertexCount, vertexStride,
                        res.iAddr, indexCount
                    );
                }

                g_ModelCache[path] = res;
                return res;
            };

            // =======================================================
            // ★ ECS 装配阶段 (加入房间和玛莎拉蒂)
            // =======================================================
            Registry registry;
            CameraControlSystem cameraControlSystem;
            CameraSystem cameraSystem;

            // 1. 摄像机
            Entity cameraEntity = registry.create();
            {
                auto& t = registry.emplace<TransformComponent>(cameraEntity);
                t.setPosition(Vector3(0.0f, 2.0f, -8.0f)); 
                
                auto& cam = registry.emplace<CameraComponent>(cameraEntity);
                cam.setFov(45.0f); cam.setAspect((float)WIDTH / (float)HEIGHT);
                cam.setzNear(0.001f); cam.setzFar(1000.0f); cam.setMain(true);

                auto& ctrl = registry.emplace<CameraControlComponent>(cameraEntity);
                ctrl.setMoveSpeed(3.0f); ctrl.setSensitivityX(0.1f); ctrl.setSensitivityY(0.1f);
                ctrl.setYaw(0.0f);
            }

            // 2. 房间/地板
            auto roomEntity = registry.create();
            auto& roomTrans = registry.emplace<TransformComponent>(roomEntity);
            roomTrans.setPosition(Vector3(0.0f, -2.0f, 0.0f));
            roomTrans.setScale(Vector3(1.0f, 1.0f, 1.0f)); // 假设你的组件有 setScale 方法，如果没有请改成直接赋值，例如 roomTrans.scale = ...
            
            registry.emplace<VulkanModelComponent>(roomEntity).setVulkanModelPath("C:/Lizeral Engine/LizeralEngine0.0.1/asset/box_with_uv.glb"); 

            // 3. 跑车阵列！
            for(int i = -1; i <= 1; i++) {
                auto carEntity = registry.create();
                auto& carTrans = registry.emplace<TransformComponent>(carEntity);
                carTrans.setPosition(Vector3(i * 5.0f, 0.0f, 0.0f));
                carTrans.setScale(Vector3(10.0f, 10.0f, 10.0f));
                
                registry.emplace<VulkanModelComponent>(carEntity).setVulkanModelPath("C:/Lizeral Engine/LizeralEngine0.0.1/asset/maserati.glb"); 
            }

            // =======================================================
            // ★ 预热机制：提前遍历 ECS 并加载模型资源
            // =======================================================
            auto view = registry.view<TransformComponent, VulkanModelComponent>();
            for (auto entity : view) {
                auto& modelComp = view.get<VulkanModelComponent>(entity);
                if (!modelComp.IsLoaded() && !modelComp.getVulkanModelPath().empty()) {
                    getOrLoadModel(modelComp.getVulkanModelPath());
                    modelComp.SetLoaded(true);
                }
            }

            // =======================================================
            // ★ 管线准备阶段 (Bindless & G-Buffer)
            // =======================================================
            std::cout << "[Sandbox] Building Global Bindless Texture Pool (" << g_GlobalImageInfos.size() << " textures)..." << std::endl;
            Lizeral::VulkanDescriptorBuilder bindlessBuilder;
            VkDescriptorSetLayout descriptorSetLayout;
            VkDescriptorPool descriptorPool;
            VkDescriptorSet descriptorSet;
            bindlessBuilder.BindImageArray(0, g_GlobalImageInfos.data(), 1024, static_cast<uint32_t>(g_GlobalImageInfos.size()), 
                                           VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT, true)
                           .Build(&vulkanDevice, descriptorSetLayout, descriptorPool, descriptorSet);

            auto createAttachment = [&](VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect) -> GBufferAttachment {
                GBufferAttachment attachment;
                attachment.format = format;

                VkImageCreateInfo imageInfo{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
                imageInfo.imageType = VK_IMAGE_TYPE_2D;
                imageInfo.extent.width = WIDTH;   
                imageInfo.extent.height = HEIGHT; 
                imageInfo.extent.depth = 1;
                imageInfo.mipLevels = 1;          
                imageInfo.arrayLayers = 1;
                imageInfo.format = format;
                imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
                imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
                imageInfo.usage = usage | VK_IMAGE_USAGE_SAMPLED_BIT; 
                imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
                imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

                VmaAllocationCreateInfo allocInfo{};
                allocInfo.usage = VMA_MEMORY_USAGE_GPU_ONLY;

                if (vmaCreateImage(vulkanDevice.GetAllocator(), &imageInfo, &allocInfo, &attachment.image, &attachment.allocation, nullptr) != VK_SUCCESS)
                    throw std::runtime_error("Failed to allocate G-Buffer image!");

                VkImageViewCreateInfo viewInfo{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
                viewInfo.image = attachment.image;
                viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
                viewInfo.format = format;
                viewInfo.subresourceRange.aspectMask = aspect;
                viewInfo.subresourceRange.baseMipLevel = 0;
                viewInfo.subresourceRange.levelCount = 1;
                viewInfo.subresourceRange.baseArrayLayer = 0;
                viewInfo.subresourceRange.layerCount = 1;

                vkCreateImageView(vulkanDevice.GetNativeDevice(), &viewInfo, nullptr, &attachment.view);
                return attachment;
            };

            GBufferAttachment gAlbedoMetallic = createAttachment(VK_FORMAT_R8G8B8A8_UNORM, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            GBufferAttachment gNormalRoughness = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            GBufferAttachment gDepth = createAttachment(VK_FORMAT_D32_SFLOAT, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT, VK_IMAGE_ASPECT_DEPTH_BIT);
            GBufferAttachment gVelocity = createAttachment(VK_FORMAT_R16G16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            GBufferAttachment gSceneColor = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

            VkSamplerCreateInfo gSamplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
            gSamplerInfo.magFilter = VK_FILTER_NEAREST; gSamplerInfo.minFilter = VK_FILTER_NEAREST;
            gSamplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
            gSamplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
            gSamplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
            gSamplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
            gSamplerInfo.maxAnisotropy = 1.0f;
            VkSampler gBufferSampler;
            vkCreateSampler(vulkanDevice.GetNativeDevice(), &gSamplerInfo, nullptr, &gBufferSampler);

            VkDescriptorImageInfo gInfos[4] = {};
            gInfos[0] = { gBufferSampler, gAlbedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };
            gInfos[1] = { gBufferSampler, gNormalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };
            gInfos[2] = { gBufferSampler, gDepth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };
            gInfos[3] = { gBufferSampler, gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };

            GBufferAttachment gHistory[2]; // 乒乓历史缓冲
            gHistory[0] = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            gHistory[1] = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

            std::cout << "[VMA Debug] SceneColor Alloc Created: " << gSceneColor.allocation << std::endl;
            std::cout << "[VMA Debug] History0 Alloc Created: " << gHistory[0].allocation << std::endl;
            std::cout << "[VMA Debug] History1 Alloc Created: " << gHistory[1].allocation << std::endl;

            VkDescriptorImageInfo taaInfos0[3] = {
                { gBufferSampler, gSceneColor.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gHistory[1].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL } // 读 1
            };
            // TAA 描述符 1 (写 History 1，读 History 0)
            VkDescriptorImageInfo taaInfos1[3] = {
                { gBufferSampler, gSceneColor.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gHistory[0].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL } // 读 0
            };

            // Blit 描述符也需要两个，跟着 Ping-Pong 走
            VkDescriptorImageInfo blitInfos0[1] = { { gBufferSampler, gHistory[0].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL } };
            VkDescriptorImageInfo blitInfos1[1] = { { gBufferSampler, gHistory[1].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL } };

            VkDescriptorSetLayout lightingDescriptorSetLayout;
            VkDescriptorPool lightDescriptorPool;
            VkDescriptorSet lightingDescriptorSet;
            Lizeral::VulkanDescriptorBuilder gBufferBuilder;
            gBufferBuilder.BindImage(0, &gInfos[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                          .BindImage(1, &gInfos[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                          .BindImage(2, &gInfos[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                          .BindImage(3, &gInfos[3], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                          .Build(&vulkanDevice, lightingDescriptorSetLayout, lightDescriptorPool, lightingDescriptorSet);

            // // =======================================================
            // // ★ 1. 构建 TAA 的两套描述符 (Set 0 和 Set 1)
            // // =======================================================
            VkDescriptorSetLayout taaDescriptorSetLayout0, taaDescriptorSetLayout1;
            VkDescriptorPool taaDescriptorPool0, taaDescriptorPool1;
            VkDescriptorSet taaSet[2];

            Lizeral::VulkanDescriptorBuilder taaBuilder0;
            taaBuilder0.BindImage(0, &taaInfos0[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                    .BindImage(1, &taaInfos0[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                    .BindImage(2, &taaInfos0[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                    .Build(&vulkanDevice, taaDescriptorSetLayout0, taaDescriptorPool0, taaSet[0]);

            Lizeral::VulkanDescriptorBuilder taaBuilder1;
            taaBuilder1.BindImage(0, &taaInfos1[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                    .BindImage(1, &taaInfos1[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                    .BindImage(2, &taaInfos1[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                    .Build(&vulkanDevice, taaDescriptorSetLayout1, taaDescriptorPool1, taaSet[1]);


            // =======================================================
            // ★ 2. 构建 Final Blit 的两套描述符 (Set 0 和 Set 1)
            // =======================================================
            VkDescriptorSetLayout blitDescriptorSetLayout0, blitDescriptorSetLayout1;
            VkDescriptorPool blitDescriptorPool0, blitDescriptorPool1;
            VkDescriptorSet blitSet[2];

            Lizeral::VulkanDescriptorBuilder blitBuilder0;
            blitBuilder0.BindImage(0, &blitInfos0[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .Build(&vulkanDevice, blitDescriptorSetLayout0, blitDescriptorPool0, blitSet[0]);

            Lizeral::VulkanDescriptorBuilder blitBuilder1;
            blitBuilder1.BindImage(0, &blitInfos1[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .Build(&vulkanDevice, blitDescriptorSetLayout1, blitDescriptorPool1, blitSet[1]);

            // 构建 Graphics 管线
            VkShaderModule meshShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_mesh.spv"));
            VkShaderModule fragShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_frag.spv"));
            VkShaderModule taskShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_task.spv"));

            VkPushConstantRange pushConstantRange{};
            pushConstantRange.stageFlags = VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT;
            pushConstantRange.offset = 0; pushConstantRange.size = sizeof(PushConstants); 

            VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
            pipelineLayoutInfo.pushConstantRangeCount = 1; pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;
            pipelineLayoutInfo.setLayoutCount = 1; pipelineLayoutInfo.pSetLayouts = &descriptorSetLayout;
            VkPipelineLayout pipelineLayout;
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &pipelineLayoutInfo, nullptr, &pipelineLayout);

            VulkanPipelineBuilder builderPipeline;
            VkPipeline graphicsPipeline = builderPipeline
                .AddShaderStage(VK_SHADER_STAGE_TASK_BIT_EXT, taskShaderModule)
                .AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShaderModule) 
                .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragShaderModule)
                .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
                .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, VK_FRONT_FACE_COUNTER_CLOCKWISE)
                .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
                .SetDepthStencil(true, true, VK_COMPARE_OP_GREATER_OR_EQUAL)
                .DisableColorBlendAttachments(3)
                .SetPipelineLayout(pipelineLayout)
                .Build(&vulkanDevice, { VK_FORMAT_R8G8B8A8_UNORM, VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16_SFLOAT }, VK_FORMAT_D32_SFLOAT); 

            // 构建 Lighting 管线
            VkShaderModule lightVertModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "lighting_vert.spv"));
            VkShaderModule lightFragModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "lighting_frag.spv"));

            VkPushConstantRange lightPushRange{};
            lightPushRange.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT; lightPushRange.offset = 0; lightPushRange.size = sizeof(LightingPushConstants);

            VkPipelineLayoutCreateInfo lightPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
            lightPipeLayoutInfo.setLayoutCount = 1; lightPipeLayoutInfo.pSetLayouts = &lightingDescriptorSetLayout;
            lightPipeLayoutInfo.pushConstantRangeCount = 1; lightPipeLayoutInfo.pPushConstantRanges = &lightPushRange; 
            VkPipelineLayout lightingPipelineLayout;
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &lightPipeLayoutInfo, nullptr, &lightingPipelineLayout);

            VulkanPipelineBuilder builderLighting;
            VkPipeline lightingPipeline = builderLighting
                .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, lightVertModule)
                .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, lightFragModule)
                .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
                .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE) 
                .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
                .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS) 
                .AddColorBlendAttachment(false) 
                .SetPipelineLayout(lightingPipelineLayout)
                .Build(&vulkanDevice, { VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED); 

            VkShaderModule taaFragModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "taa_frag.spv"));

            VkPipelineLayoutCreateInfo taaPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
            taaPipeLayoutInfo.setLayoutCount = 1; 
            // 提示：就用 0 号 Layout，因为 0 和 1 结构完全一样
            taaPipeLayoutInfo.pSetLayouts = &taaDescriptorSetLayout0; 
            VkPipelineLayout taaPipelineLayout;
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &taaPipeLayoutInfo, nullptr, &taaPipelineLayout);

            VulkanPipelineBuilder builderTaa;
            VkPipeline taaPipeline = builderTaa
                .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, lightVertModule) // ★ 复用全屏三角形
                .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, taaFragModule)
                .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
                .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE) 
                .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
                .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)         // 后处理不需要深度测试
                .AddColorBlendAttachment(false) 
                .SetPipelineLayout(taaPipelineLayout)
                // ★ TAA 输出到 History 缓冲，格式必须一致；不需要 Depth 附件
                .Build(&vulkanDevice, { VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED); 


            // =======================================================
            // ★ 4. 构建 Final Blit 管线 (上屏 Pipeline)
            // =======================================================
            VkShaderModule blitFragModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "blit_frag.spv"));

            VkPipelineLayoutCreateInfo blitPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
            blitPipeLayoutInfo.setLayoutCount = 1; 
            blitPipeLayoutInfo.pSetLayouts = &blitDescriptorSetLayout0; 
            VkPipelineLayout blitPipelineLayout;
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &blitPipeLayoutInfo, nullptr, &blitPipelineLayout);

            VulkanPipelineBuilder builderBlit;
            VkPipeline blitPipeline = builderBlit
                .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, lightVertModule) // ★ 再次复用
                .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, blitFragModule)
                .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
                .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE) 
                .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
                .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS) 
                .AddColorBlendAttachment(false) 
                .SetPipelineLayout(blitPipelineLayout)
                .Build(&vulkanDevice, { renderer.GetSwapchainFormat() }, VK_FORMAT_D32_SFLOAT);

            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), lightVertModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), lightFragModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), fragShaderModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), meshShaderModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), taskShaderModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), taaFragModule,   nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), blitFragModule,   nullptr);

            std::cout << "\n[Sandbox] Dynamic Rendering & ECS Ready! Hold RMB + WASD to move." << std::endl;

            auto CmdDrawMeshTasksEXT = (PFN_vkCmdDrawMeshTasksEXT)vkGetDeviceProcAddr(vulkanDevice.GetNativeDevice(), "vkCmdDrawMeshTasksEXT");
            auto vkCmdBeginRendering = (PFN_vkCmdBeginRendering)vkGetDeviceProcAddr(vulkanDevice.GetNativeDevice(), "vkCmdBeginRendering");
            auto vkCmdEndRendering = (PFN_vkCmdEndRendering)vkGetDeviceProcAddr(vulkanDevice.GetNativeDevice(), "vkCmdEndRendering");

            bool firstFrame = true;

            Matrix4x4 prevViewProj;
            std::unordered_map<Entity, Matrix4x4> prevModelMats;
            bool isFirstFrameRun = true;

            static uint32_t frameIndex = 0;

            float lastTime = glfwGetTime();

            std::cout << "[Sandbox] Building TLAS..." << std::endl;

            std::vector<VkAccelerationStructureInstanceKHR> tlasInstances;
            uint32_t customInstanceId = 0;

            auto renderView = registry.view<TransformComponent, VulkanModelComponent>();
            for (auto entity : renderView) {
                auto& transform = renderView.get<TransformComponent>(entity);
                auto& modelComp = renderView.get<VulkanModelComponent>(entity);
                if (!modelComp.IsLoaded()) continue;

                const Lizeral::VulkanModelResource& res = g_ModelCache[modelComp.getVulkanModelPath()];
                if (!res.blas) continue;

                Matrix4x4 modelMat = transform.getMatrix();
                VkTransformMatrixKHR vkTransform{};
                vkTransform.matrix[0][0] = modelMat[0][0]; vkTransform.matrix[0][1] = modelMat[0][1]; vkTransform.matrix[0][2] = modelMat[0][2]; vkTransform.matrix[0][3] = modelMat[0][3];
                vkTransform.matrix[1][0] = modelMat[1][0]; vkTransform.matrix[1][1] = modelMat[1][1]; vkTransform.matrix[1][2] = modelMat[1][2]; vkTransform.matrix[1][3] = modelMat[1][3];
                vkTransform.matrix[2][0] = modelMat[2][0]; vkTransform.matrix[2][1] = modelMat[2][1]; vkTransform.matrix[2][2] = modelMat[2][2]; vkTransform.matrix[2][3] = modelMat[2][3];

                VkAccelerationStructureInstanceKHR instance{};
                instance.transform = vkTransform;
                instance.instanceCustomIndex = customInstanceId++; 
                instance.mask = 0xFF; 
                instance.instanceShaderBindingTableRecordOffset = 0; 
                instance.flags = VK_GEOMETRY_INSTANCE_TRIANGLE_FACING_CULL_DISABLE_BIT_KHR; 
                instance.accelerationStructureReference = res.blas->GetDeviceAddress();

                tlasInstances.push_back(instance);
            }

            if (!tlasInstances.empty()) {
                // 借用你写好的 VulkanCommandBuffer 类来做一次性提交
                Lizeral::VulkanCommandBuffer tlasCmd(&vulkanDevice, &resourceCommandPool);
                tlasCmd.Begin(VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT);

                tlas.Build(tlasCmd.GetNativeBuffer(), tlasInstances);

                // 加一个内存屏障，确保 TLAS 完全建好
                VkMemoryBarrier memoryBarrier{};
                memoryBarrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
                memoryBarrier.srcAccessMask = VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR;
                memoryBarrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
                vkCmdPipelineBarrier(tlasCmd.GetNativeBuffer(), 
                    VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR, 
                    VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 
                    0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);

                tlasCmd.End();
                tlasCmd.SubmitAndIdle(); // 暴力等待 GPU 彻底建完 TLAS
                std::cout << "[Sandbox] TLAS Built Successfully!" << std::endl;
            }

            std::cout<<"enter main loop"<<std::endl;

            while (!glfwWindowShouldClose(window)) {
                
                float currentTime = glfwGetTime();
                float dt = currentTime - lastTime;
                lastTime = currentTime;
                frameIndex++;

                uint32_t jitterIndex = frameIndex % 8 + 1; 
                float jitterX = CreateHaltonSequence(jitterIndex, 2) - 0.5f;
                float jitterY = CreateHaltonSequence(jitterIndex, 3) - 0.5f;

                // 2. 将亚像素偏移转换到 NDC 空间 [-1, 1]
                // 乘以 2 是因为 NDC 宽度是 2，除以屏幕宽高是因为抖动幅度刚好是一个像素
                float ndcJitterX = (jitterX * 2.0f) / static_cast<float>(WIDTH);
                float ndcJitterY = (jitterY * 2.0f) / static_cast<float>(HEIGHT);

                Input::GetInstance().Tick(); 
                if (Input::GetInstance().GetKey(Key::ESC)) glfwSetWindowShouldClose(window, true);
                glfwPollEvents();

                cameraControlSystem.Tick(dt, registry);
                cameraSystem.Tick(registry); 

                Matrix4x4 viewMat, projMatUnjittered, projMatJittered;
                auto cameraView = registry.view<TransformComponent, CameraComponent>();
                for (auto entity : cameraView) {
                    auto& camera = cameraView.get<CameraComponent>(entity);
                    camera.setProjectionMatrix(camera.BuildPerspectiveInfiniteReverseZ(camera.getFov(),camera.getAspect(),camera.getzNear()));
                    viewMat = camera.getViewMatrix(); 
                    // projMat = camera.BuildPerspectiveInfiniteReverseZ(camera.getFov(),camera.getAspect(),camera.getzNear());

                    projMatUnjittered = camera.BuildPerspectiveInfiniteReverseZ(camera.getFov(), camera.getAspect(), camera.getzNear());
                    projMatUnjittered[1][1] *= -1.0f; // Vulkan Y 轴补丁
                    
                    projMatJittered = projMatUnjittered;

                    projMatJittered[0][2] += ndcJitterX;
                    projMatJittered[1][2] += ndcJitterY;
                    break; 
                }

                // 在进入 Lighting Pass 前，计算乒乓索引
                uint32_t ping = frameIndex % 2;       // 当前帧，用于写入
                uint32_t pong = (frameIndex + 1) % 2; // 上一帧，用于读取

                VkCommandBuffer cmd = renderer.BeginFrame();

                if (cmd != VK_NULL_HANDLE) {
                    
                    // 状态转换：准备写入 G-Buffer
                    VkImageLayout currentLayout = firstFrame ? VK_IMAGE_LAYOUT_UNDEFINED : VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
                    transitionImageLayout(cmd, gAlbedoMetallic.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gNormalRoughness.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gDepth.image, currentLayout, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);
                    transitionImageLayout(cmd, gVelocity.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    firstFrame = false;

                    // 配置 G-Buffer
                    VkRenderingAttachmentInfo colorAttachments[3] = {};
                    colorAttachments[0].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    colorAttachments[0].imageView = gAlbedoMetallic.view;
                    colorAttachments[0].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                    colorAttachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; colorAttachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    colorAttachments[0].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};

                    colorAttachments[1].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    colorAttachments[1].imageView = gNormalRoughness.view;
                    colorAttachments[1].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                    colorAttachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; colorAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    colorAttachments[1].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};

                    colorAttachments[2].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    colorAttachments[2].imageView = gVelocity.view;
                    colorAttachments[2].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                    colorAttachments[2].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; 
                    colorAttachments[2].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    colorAttachments[2].clearValue.color = {0.0f, 0.0f, 0.0f, 0.0f};

                    VkRenderingAttachmentInfo depthAttachment{};
                    depthAttachment.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    depthAttachment.imageView = gDepth.view;
                    depthAttachment.imageLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
                    depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    depthAttachment.clearValue.depthStencil = {0.0f, 0};

                    VkRenderingInfo renderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO};
                    renderInfo.renderArea = {{0, 0}, {WIDTH, HEIGHT}};
                    renderInfo.layerCount = 1; renderInfo.colorAttachmentCount = 3;
                    renderInfo.pColorAttachments = colorAttachments; renderInfo.pDepthAttachment = &depthAttachment;

                    vkCmdBeginRendering(cmd, &renderInfo);

                    VkViewport viewport{};
                    viewport.width = static_cast<float>(WIDTH); viewport.height = static_cast<float>(HEIGHT); viewport.maxDepth = 1.0f;
                    vkCmdSetViewport(cmd, 0, 1, &viewport);

                    VkRect2D scissor{};
                    scissor.extent = {static_cast<uint32_t>(WIDTH), static_cast<uint32_t>(HEIGHT)};
                    vkCmdSetScissor(cmd, 0, 1, &scissor);
                    
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, graphicsPipeline);
                    // 绑定全局宇宙贴图数组！
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout, 0, 1, &descriptorSet, 0, nullptr);

                    Matrix4x4 currentViewProjJittered = projMatJittered * viewMat;
                    Matrix4x4 currentViewProjUnjittered = projMatUnjittered * viewMat;

                    if (isFirstFrameRun) prevViewProj = currentViewProjUnjittered;
                    
                    // ====================================================================
                    // ★ 数据驱动渲染发车
                    // ====================================================================
                    for (auto entity : renderView) {
                        auto& transform = renderView.get<TransformComponent>(entity);
                        auto& modelComp = renderView.get<VulkanModelComponent>(entity);

                        if (!modelComp.IsLoaded()) continue;

                        const Lizeral::VulkanModelResource& res = g_ModelCache[modelComp.getVulkanModelPath()];

                        Matrix4x4 currentModel = transform.getMatrix();
                        Matrix4x4 prevModel = isFirstFrameRun ? currentModel : prevModelMats[entity];

                        PushConstants pushData{};
                        // pushData.mvp = (currentViewProjJittered * currentModel).transpose();
                        pushData.mvp = (currentViewProjUnjittered * currentModel).transpose();
                        pushData.model = transform.getMatrix().transpose(); 
                        pushData.prevMvp = (prevViewProj * prevModel).transpose();
                        pushData.vertexBuffer = res.vAddr;
                        pushData.meshletBuffer = res.mAddr;
                        pushData.indexBuffer = res.iAddr;
                        pushData.boundsBuffer = res.bAddr; 
                        pushData.materialBuffer = res.matAddr;
                        pushData.totalMeshlets = res.totalMeshlets; 
                        pushData.textureOffset = res.textureOffset;
                        pushData.jitter = Vector2(ndcJitterX, ndcJitterY);

                        vkCmdPushConstants(cmd, pipelineLayout, VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(PushConstants), &pushData);
                        
                        prevModelMats[entity] = currentModel;

                        uint32_t taskGroupCount = (res.totalMeshlets + 63) / 64;
                        CmdDrawMeshTasksEXT(cmd, taskGroupCount, 1, 1);
                    }

                    vkCmdEndRendering(cmd);

                    prevViewProj = currentViewProjUnjittered;
                    isFirstFrameRun = false;

                    // 计算乒乓缓冲索引：这帧写 ping，读 pong
                    uint32_t ping = frameIndex % 2;       
                    uint32_t pong = (frameIndex + 1) % 2; 


                    // ★ 初始化历史贴图的起始状态 (仅在程序运行的第一帧执行)
                    static bool isFirstTAAFrame = true;
                    if (isFirstTAAFrame) {
                        transitionImageLayout(cmd, gSceneColor.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                        transitionImageLayout(cmd, gHistory[0].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                        transitionImageLayout(cmd, gHistory[1].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                        isFirstTAAFrame = false;
                    }
                    

                    // 状态转换：准备被光照 Shader 读取
                    transitionImageLayout(cmd, gAlbedoMetallic.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gNormalRoughness.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gDepth.image, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);
                    transitionImageLayout(cmd, gVelocity.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

                    // ====================================================================
                    // ★ Pass 2: 全屏光照计算通道 (改道：画到离屏贴图 gSceneColor)
                    // ====================================================================
                    // 将 SceneColor 从只读转换为可写
                    transitionImageLayout(cmd, gSceneColor.image, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

                    VkRenderingAttachmentInfo sceneColorAttachment{};
                    sceneColorAttachment.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    sceneColorAttachment.imageView = gSceneColor.view;
                    sceneColorAttachment.imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                    sceneColorAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; 
                    sceneColorAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    sceneColorAttachment.clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};

                    VkRenderingInfo lightRenderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO};
                    lightRenderInfo.renderArea = {{0, 0}, {WIDTH, HEIGHT}};
                    lightRenderInfo.layerCount = 1; 
                    lightRenderInfo.colorAttachmentCount = 1;
                    lightRenderInfo.pColorAttachments = &sceneColorAttachment;

                    vkCmdBeginRendering(cmd, &lightRenderInfo); 
                    
                    vkCmdSetViewport(cmd, 0, 1, &viewport);
                    vkCmdSetScissor(cmd, 0, 1, &scissor);

                    LightingPushConstants lightPc{};
                    Matrix4x4 vpJittered = projMatJittered * viewMat;
                    lightPc.invViewProj = vpJittered.inverse().transpose(); 
                    lightPc.viewProj = vpJittered.transpose();
                    lightPc.cameraPos = registry.get<TransformComponent>(cameraEntity).getPosition(); 

                    vkCmdPushConstants(cmd, lightingPipelineLayout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(LightingPushConstants), &lightPc);
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, lightingPipeline);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, lightingPipelineLayout, 0, 1, &lightingDescriptorSet, 0, nullptr);
                    
                    vkCmdDraw(cmd, 3, 1, 0, 0);

                    vkCmdEndRendering(cmd);

                    // // ====================================================================
                    // // ★ Pass 3: TAA 时空融合通道 (写到 History[ping])
                    // // ====================================================================
                    // // SceneColor 刚画完，转回只读给 TAA 用
                    transitionImageLayout(cmd, gSceneColor.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    // 准备把 TAA 结果写入到当前的 History[ping] 中，所以把它变成可写
                    transitionImageLayout(cmd, gHistory[ping].image, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

                    VkRenderingAttachmentInfo taaAttachment{};
                    taaAttachment.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    taaAttachment.imageView = gHistory[ping].view;
                    taaAttachment.imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                    taaAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE; // 画面会全覆盖，所以不用 Clear，节约性能
                    taaAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;

                    VkRenderingInfo taaRenderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO};
                    taaRenderInfo.renderArea = {{0, 0}, {WIDTH, HEIGHT}};
                    taaRenderInfo.layerCount = 1; 
                    taaRenderInfo.colorAttachmentCount = 1;
                    taaRenderInfo.pColorAttachments = &taaAttachment;

                    vkCmdBeginRendering(cmd, &taaRenderInfo);
                    vkCmdSetViewport(cmd, 0, 1, &viewport);
                    vkCmdSetScissor(cmd, 0, 1, &scissor);

                    // 绑定你在上一步创建的 TAA 管线和动态 Ping-Pong 描述符
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, taaPipeline);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, taaPipelineLayout, 0, 1, &taaSet[ping], 0, nullptr);
                    vkCmdDraw(cmd, 3, 1, 0, 0);
                    vkCmdEndRendering(cmd);

                    // ====================================================================
                    // ★ Pass 4: Final Blit 通道 (画到 Swapchain 显示器上)
                    // ====================================================================
                    // 刚刚融合好的 History[ping] 转换成只读，准备被读取上屏
                    transitionImageLayout(cmd, gHistory[ping].image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

                    // 这里使用你原有的 renderer.BeginRendering，它会自动绑定 Swapchain 的图像作为目标！
                    renderer.BeginRendering(cmd); 
                    
                    vkCmdSetViewport(cmd, 0, 1, &viewport);
                    vkCmdSetScissor(cmd, 0, 1, &scissor);

                    // 绑定 Blit 搬运管线
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, blitPipeline);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, blitPipelineLayout, 0, 1, &blitSet[ping], 0, nullptr);
                    
                    vkCmdDraw(cmd, 3, 1, 0, 0);

                    renderer.EndRendering(cmd);
                    renderer.EndFrame();
                }
            }

            std::cout << "--- EXITING PROGRAM ---" << std::endl;
            
            vkDeviceWaitIdle(vulkanDevice.GetNativeDevice());

            // ====================================================================
            // ★ 清理与善后
            // ====================================================================
            std::cout<<"[Pipeline] start destroying Pipeline"<<std::endl;
            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), blitPipeline, nullptr);
            vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), blitPipelineLayout, nullptr);

            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), taaPipeline, nullptr);
            vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), taaPipelineLayout, nullptr);

            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), lightingPipeline, nullptr);
            vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), lightingPipelineLayout, nullptr);

            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), graphicsPipeline, nullptr);
            vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), pipelineLayout, nullptr);

            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), descriptorPool, nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), lightDescriptorPool, nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), taaDescriptorPool0, nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), taaDescriptorPool1, nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), blitDescriptorPool0, nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), blitDescriptorPool1, nullptr);

            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), descriptorSetLayout, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), lightingDescriptorSetLayout, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), blitDescriptorSetLayout0, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), blitDescriptorSetLayout1, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), taaDescriptorSetLayout0, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), taaDescriptorSetLayout1, nullptr);

            vkDestroySampler(vulkanDevice.GetNativeDevice(), gBufferSampler, nullptr);

            std::cout<< "[VMA Debug] start clear BDA buffer"<<std::endl;
            
            // 清理所有的 BDA 缓冲
            for (auto& alloc : g_AllocatedBuffers) {
                vkDestroyBuffer(vulkanDevice.GetNativeDevice(), alloc.buffer, nullptr);
                vkFreeMemory(vulkanDevice.GetNativeDevice(), alloc.memory, nullptr);
            }


            // 清理 G-Buffer
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gAlbedoMetallic.view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gAlbedoMetallic.image, gAlbedoMetallic.allocation);

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gNormalRoughness.view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gNormalRoughness.image, gNormalRoughness.allocation);

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gDepth.view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gDepth.image, gDepth.allocation);

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gVelocity.view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gVelocity.image, gVelocity.allocation);

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gSceneColor.view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gSceneColor.image, gSceneColor.allocation);

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gHistory[0].view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gHistory[0].image, gHistory[0].allocation);

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gHistory[1].view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gHistory[1].image, gHistory[1].allocation);
            // 全局贴图池 g_GlobalTextures 会因为 unique_ptr 自动析构
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