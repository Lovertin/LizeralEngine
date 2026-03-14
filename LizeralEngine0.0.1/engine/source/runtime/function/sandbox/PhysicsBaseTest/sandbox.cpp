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
#include "runtime/function/ecs/components/Model/VulkanModelComponent.h" 
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

struct RTInstanceDesc {
    uint64_t vertexBuffer;
    uint64_t indexBuffer;
    uint64_t materialBuffer;
    uint32_t textureOffset;
    uint32_t padding[3]; // 保持 16 字节对齐
};

struct LightingPushConstants {
    Lizeral::Matrix4x4 invViewProj; 
    Lizeral::Matrix4x4 viewProj;
    Lizeral::Vector3 cameraPos;   
    float padding;  
    uint32_t frameIndex;
    uint32_t padding2;           // 强制对齐
    uint64_t instanceDescAddr;   // 新增：传递全局台账的 BDA 指针  
    Vector3 lightDir;   
    float lightPadding;  
    Vector3 lightColor;
    float lightIntensity;             
};

struct GBufferAttachment {
    VkImage image = VK_NULL_HANDLE;
    VmaAllocation allocation = VK_NULL_HANDLE;
    VkImageView view = VK_NULL_HANDLE;
    VkFormat format;
};

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
    } else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
        barrier.srcAccessMask = 0;
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
        destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    } else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL) {
        barrier.srcAccessMask = 0; 
        barrier.dstAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
        sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT; 
        destinationStage = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT;
    } else if (oldLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
        barrier.srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT; 
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        sourceStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT; 
        destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    } else if (oldLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
        barrier.srcAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT; 
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        sourceStage = VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT; 
        destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    } else if (oldLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL) {
        barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT; 
        barrier.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
        sourceStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT; 
        destinationStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    } else if (oldLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL) {
        barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT;
        barrier.dstAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
        sourceStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
        destinationStage = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT; 
    } else {
        std::cerr << "Unsupported transition: " << oldLayout << " -> " << newLayout << std::endl;
        throw std::invalid_argument("Unsupported layout transition!");
    }

    vkCmdPipelineBarrier(cmd, sourceStage, destinationStage, 0, 0, nullptr, 0, nullptr, 1, &barrier);
}

float CreateHaltonSequence(uint32_t index, uint32_t base) {
    float f = 1.0f; float result = 0.0f; uint32_t i = index;
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
        window = glfwCreateWindow(WIDTH, HEIGHT, "Lizeral Engine - Hardware Ray Tracing", nullptr, nullptr);

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

            std::unordered_map<std::string, Lizeral::VulkanModelResource> g_ModelCache;
            std::vector<std::unique_ptr<Lizeral::VulkanTexture>> g_GlobalTextures; 
            std::vector<VkDescriptorImageInfo> g_GlobalImageInfos;
            std::vector<BufferAllocation> g_AllocatedBuffers;

            VulkanTLAS tlas(&vulkanDevice, 2);

            auto getOrLoadModel = [&](const std::string& path) -> Lizeral::VulkanModelResource {
                if (g_ModelCache.find(path) != g_ModelCache.end()) return g_ModelCache[path];
                std::cout << "[RenderSystem] Loading new model to GPU: " << path << std::endl;
                
                uint32_t currentTexOffset = static_cast<uint32_t>(g_GlobalTextures.size());
                MeshletModelBuilder builder;
                if (!builder.LoadAndSliceModel(path, currentTexOffset)) throw std::runtime_error("Failed to load GLB!");
                
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
                
                g_AllocatedBuffers.push_back({vBuf, vMem}); g_AllocatedBuffers.push_back({mBuf, mMem});
                g_AllocatedBuffers.push_back({iBuf, iMem}); g_AllocatedBuffers.push_back({bBuf, bMem});
                g_AllocatedBuffers.push_back({matBuf, matMem});

                res.totalMeshlets = builder.GetMeshlets().size();
                res.textureOffset = currentTexOffset;
                res.textureCount = builder.GetAllTextures().size();

                const auto& vertices = builder.GetVertices();
                const auto& microIndices = builder.GetMicroIndices();
                const auto& meshlets = builder.GetMeshlets();

                // 1. 翻译局部索引为全局索引
                std::vector<uint32_t> globalIndices;
                globalIndices.reserve(microIndices.size());

                for (const auto& m : meshlets) {
                    for (uint32_t i = 0; i < m.triangleCount * 3; i++) {
                        // 真正的全局顶点索引 = Meshlet的顶点基础偏移 + 内部微索引
                        globalIndices.push_back(m.vertexOffset + microIndices[m.triangleOffset + i]); 
                    }
                }

                uint32_t vertexCount = static_cast<uint32_t>(vertices.size());
                uint32_t indexCount = static_cast<uint32_t>(globalIndices.size());
                uint32_t vertexStride = vertices.empty() ? 0 : static_cast<uint32_t>(sizeof(vertices[0]));

                if (vertexCount > 0 && indexCount > 0) {
                    // 2. 为光追专门上传一份全局索引缓冲 (BLAS 专用！)
                    VkBuffer blasIBuf;
                    VkDeviceMemory blasIMem;
                    uint64_t blasIAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), globalIndices, blasIBuf, blasIMem);
                    res.globalIAddr = blasIAddr;
                    
                    // 登记到销毁清单
                    g_AllocatedBuffers.push_back({blasIBuf, blasIMem});

                    std::cout << "[VulkanBLAS] Triggering BLAS build for: " << path << " with Global Indices" << std::endl;
                    
                    // 3. 实例化并构建 BLAS，★注意：这里传入的是 blasIAddr，而不是 res.iAddr！
                    res.blas = std::make_shared<VulkanBLAS>(
                        &vulkanDevice, 
                        &resourceCommandPool, 
                        res.vAddr, vertexCount, vertexStride,
                        blasIAddr, indexCount // 变了！
                    );
                }

                g_ModelCache[path] = res;
                return res;
            };

            // =======================================================
            // ★ ECS 装配阶段
            // =======================================================
            Registry registry;
            CameraControlSystem cameraControlSystem;
            CameraSystem cameraSystem;

            Entity cameraEntity = registry.create();
            {
                auto& t = registry.emplace<TransformComponent>(cameraEntity);
                t.setPosition(Vector3(0.0f, 2.0f, 2.0f)); 
                auto& cam = registry.emplace<CameraComponent>(cameraEntity);
                cam.setFov(45.0f); cam.setAspect((float)WIDTH / (float)HEIGHT);
                cam.setzNear(0.001f); cam.setzFar(1000.0f); cam.setMain(true);
                auto& ctrl = registry.emplace<CameraControlComponent>(cameraEntity);
                ctrl.setMoveSpeed(3.0f); ctrl.setSensitivityX(0.1f); ctrl.setSensitivityY(0.1f);
                ctrl.setSpeedMutiplier(10.0f);
            }

            auto roomEntity = registry.create();
            auto& roomTrans = registry.emplace<TransformComponent>(roomEntity);
            roomTrans.setPosition(Vector3(0.0f, -1.0f, 10.0f));
            roomTrans.setScale(Vector3(0.05f, 0.05f, 0.05f)); 
            registry.emplace<VulkanModelComponent>(roomEntity).setVulkanModelPath("C:/Lizeral Engine/LizeralEngine0.0.1/asset/scene_without_window.glb"); 

            auto sponza = registry.create();
            auto& sponzaTrans = registry.emplace<TransformComponent>(sponza);
            sponzaTrans.setPosition(Vector3(100.0f, -1.0f, 10.0f));
            sponzaTrans.setScale(Vector3(1.0f, 1.1f, 1.1f)); 
            registry.emplace<VulkanModelComponent>(sponza).setVulkanModelPath("C:/Lizeral Engine/LizeralEngine0.0.1/asset/Sponza/Sponza/glTF/Sponza.gltf"); 

            // 预热加载
            auto view = registry.view<TransformComponent, VulkanModelComponent>();
            for (auto entity : view) {
                auto& modelComp = view.get<VulkanModelComponent>(entity);
                if (!modelComp.IsLoaded() && !modelComp.getVulkanModelPath().empty()) {
                    getOrLoadModel(modelComp.getVulkanModelPath());
                    modelComp.SetLoaded(true);
                }
            }

            std::cout << "[Sandbox] Pre-building Ping-Pong TLAS..." << std::endl;

            std::vector<VkAccelerationStructureInstanceKHR> initialTlasInstances;
            uint32_t initInstanceId = 0;

            for (auto entity : view) {
                auto& transform = view.get<TransformComponent>(entity);
                auto& modelComp = view.get<VulkanModelComponent>(entity);
                if (!modelComp.IsLoaded()) continue;
                const auto& res = g_ModelCache[modelComp.getVulkanModelPath()];
                if (!res.blas) continue;

                Matrix4x4 modelMat = transform.getMatrix();
                VkTransformMatrixKHR vkTransform{};
                vkTransform.matrix[0][0] = modelMat[0][0]; vkTransform.matrix[0][1] = modelMat[0][1]; vkTransform.matrix[0][2] = modelMat[0][2]; vkTransform.matrix[0][3] = modelMat[0][3];
                vkTransform.matrix[1][0] = modelMat[1][0]; vkTransform.matrix[1][1] = modelMat[1][1]; vkTransform.matrix[1][2] = modelMat[1][2]; vkTransform.matrix[1][3] = modelMat[1][3];
                vkTransform.matrix[2][0] = modelMat[2][0]; vkTransform.matrix[2][1] = modelMat[2][1]; vkTransform.matrix[2][2] = modelMat[2][2]; vkTransform.matrix[2][3] = modelMat[2][3];

                VkAccelerationStructureInstanceKHR instance{};
                instance.transform = vkTransform;
                instance.instanceCustomIndex = initInstanceId++; 
                instance.mask = 0xFF; 
                instance.instanceShaderBindingTableRecordOffset = 0; 
                instance.flags = VK_GEOMETRY_INSTANCE_TRIANGLE_FACING_CULL_DISABLE_BIT_KHR; 
                instance.accelerationStructureReference = res.blas->GetDeviceAddress();
                initialTlasInstances.push_back(instance);
            }

            if (!initialTlasInstances.empty()) {
                Lizeral::VulkanCommandBuffer tlasCmd(&vulkanDevice, &resourceCommandPool);
                tlasCmd.Begin(VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT);

                // 建第 0 帧
                tlas.Build(tlasCmd.GetNativeBuffer(), 0, initialTlasInstances);
                VkMemoryBarrier barrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
                barrier.srcAccessMask = VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR;
                barrier.dstAccessMask = VK_ACCESS_ACCELERATION_STRUCTURE_READ_BIT_KHR | VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR;
                vkCmdPipelineBarrier(tlasCmd.GetNativeBuffer(), VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR, VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR, 0, 1, &barrier, 0, nullptr, 0, nullptr);

                // 建第 1 帧
                tlas.Build(tlasCmd.GetNativeBuffer(), 1, initialTlasInstances);
                barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
                vkCmdPipelineBarrier(tlasCmd.GetNativeBuffer(), VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0, 1, &barrier, 0, nullptr, 0, nullptr);

                tlasCmd.End();
                tlasCmd.SubmitAndIdle(); // ★ 死等！现在 TLAS 有真实的 Handle 了！
                std::cout << "[Sandbox] Ping-Pong TLAS Pre-built Successfully!" << std::endl;
            }

            // =======================================================
            // ★ 管线准备阶段
            // =======================================================
            Lizeral::VulkanDescriptorBuilder bindlessBuilder;
            VkDescriptorSetLayout descriptorSetLayout; VkDescriptorPool descriptorPool; VkDescriptorSet descriptorSet;
            bindlessBuilder.BindImageArray(0, g_GlobalImageInfos.data(), 1024, static_cast<uint32_t>(g_GlobalImageInfos.size()), 
                                           VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT, true)
                           .Build(&vulkanDevice, descriptorSetLayout, descriptorPool, descriptorSet);

            auto createAttachment = [&](VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect) -> GBufferAttachment {
                GBufferAttachment attachment; attachment.format = format;
                VkImageCreateInfo imageInfo{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
                imageInfo.imageType = VK_IMAGE_TYPE_2D; imageInfo.extent = {WIDTH, HEIGHT, 1};
                imageInfo.mipLevels = 1; imageInfo.arrayLayers = 1; imageInfo.format = format;
                imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL; imageInfo.usage = usage | VK_IMAGE_USAGE_SAMPLED_BIT; 
                imageInfo.samples = VK_SAMPLE_COUNT_1_BIT; imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
                VmaAllocationCreateInfo allocInfo{}; allocInfo.usage = VMA_MEMORY_USAGE_GPU_ONLY;
                vmaCreateImage(vulkanDevice.GetAllocator(), &imageInfo, &allocInfo, &attachment.image, &attachment.allocation, nullptr);

                VkImageViewCreateInfo viewInfo{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
                viewInfo.image = attachment.image; viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D; viewInfo.format = format;
                viewInfo.subresourceRange.aspectMask = aspect; viewInfo.subresourceRange.levelCount = 1; viewInfo.subresourceRange.layerCount = 1;
                vkCreateImageView(vulkanDevice.GetNativeDevice(), &viewInfo, nullptr, &attachment.view);
                return attachment;
            };

            GBufferAttachment gAlbedoMetallic = createAttachment(VK_FORMAT_R8G8B8A8_UNORM, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            GBufferAttachment gNormalRoughness = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            GBufferAttachment gDepth = createAttachment(VK_FORMAT_D32_SFLOAT, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT, VK_IMAGE_ASPECT_DEPTH_BIT);
            GBufferAttachment gVelocity = createAttachment(VK_FORMAT_R16G16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            // GBufferAttachment gSceneColor = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

            GBufferAttachment gDirectLight = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            GBufferAttachment gNoisyGI = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            GBufferAttachment gDenoisedGI = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

            VkSamplerCreateInfo gSamplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
            gSamplerInfo.magFilter = VK_FILTER_LINEAR; gSamplerInfo.minFilter = VK_FILTER_LINEAR;
            gSamplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
            gSamplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE; gSamplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE; gSamplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
            VkSampler gBufferSampler;
            vkCreateSampler(vulkanDevice.GetNativeDevice(), &gSamplerInfo, nullptr, &gBufferSampler);

            VkDescriptorImageInfo gInfos[4] = {
                { gBufferSampler, gAlbedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gNormalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gDepth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
            };

            VkDescriptorImageInfo denoiseInfos[5] = {
                { gBufferSampler, gNoisyGI.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gNormalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gDepth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gDirectLight.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },   // 新增
                { gBufferSampler, gAlbedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL } // 新增
            };

            GBufferAttachment gHistory[2]; 
            gHistory[0] = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
            gHistory[1] = createAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

            // =======================================================
            // ★ 给第 0 帧用的 TAA 输入 (读上一帧的 History[1])
            // =======================================================
            VkDescriptorImageInfo taaInfos0[3] = { 
                { gBufferSampler, gDenoisedGI.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gHistory[1].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }, 
                { gBufferSampler, gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }    
            };

            // =======================================================
            // ★ 给第 1 帧用的 TAA 输入 (读上一帧的 History[0])
            // =======================================================
            VkDescriptorImageInfo taaInfos1[3] = { 
                { gBufferSampler, gDenoisedGI.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { gBufferSampler, gHistory[0].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }, 
                { gBufferSampler, gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }    
            };

            VkDescriptorImageInfo blitInfos0[1] = { { gBufferSampler, gHistory[0].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL } };

            VkDescriptorImageInfo blitInfos1[1] = { { gBufferSampler, gHistory[1].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL } };

            // =======================================================
            // ★ 光照 Descriptor 绑定 (此时 Handle 绝对不为空！)
            // =======================================================
            VkDescriptorSetLayout lightingLayouts[2]; VkDescriptorPool lightPools[2]; VkDescriptorSet lightingSets[2];
            VkAccelerationStructureKHR tlasHandle0 = tlas.GetHandle(0);
            VkAccelerationStructureKHR tlasHandle1 = tlas.GetHandle(1);

            Lizeral::VulkanDescriptorBuilder gBufferBuilder0;
            gBufferBuilder0.BindImage(0, &gInfos[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .BindImage(1, &gInfos[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .BindImage(2, &gInfos[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .BindImage(3, &gInfos[3], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .BindAccelerationStructure(4, &tlasHandle0, VK_SHADER_STAGE_FRAGMENT_BIT) 
                        .Build(&vulkanDevice, lightingLayouts[0], lightPools[0], lightingSets[0]);

            Lizeral::VulkanDescriptorBuilder gBufferBuilder1;
            gBufferBuilder1.BindImage(0, &gInfos[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .BindImage(1, &gInfos[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .BindImage(2, &gInfos[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .BindImage(3, &gInfos[3], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                        .BindAccelerationStructure(4, &tlasHandle1, VK_SHADER_STAGE_FRAGMENT_BIT) 
                        .Build(&vulkanDevice, lightingLayouts[1], lightPools[1], lightingSets[1]);

            VkDescriptorSetLayout taaLayout0, taaLayout1; VkDescriptorPool taaPool0, taaPool1; VkDescriptorSet taaSet[2];
            Lizeral::VulkanDescriptorBuilder()
                .BindImage(0, &taaInfos0[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(1, &taaInfos0[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(2, &taaInfos0[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .Build(&vulkanDevice, taaLayout0, taaPool0, taaSet[0]);

            Lizeral::VulkanDescriptorBuilder()
                .BindImage(0, &taaInfos1[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(1, &taaInfos1[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(2, &taaInfos1[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .Build(&vulkanDevice, taaLayout1, taaPool1, taaSet[1]);

            VkDescriptorSetLayout blitLayout0, blitLayout1; VkDescriptorPool blitPool0, blitPool1; VkDescriptorSet blitSet[2];
            Lizeral::VulkanDescriptorBuilder().BindImage(0, &blitInfos0[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT).Build(&vulkanDevice, blitLayout0, blitPool0, blitSet[0]);
            Lizeral::VulkanDescriptorBuilder().BindImage(0, &blitInfos1[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT).Build(&vulkanDevice, blitLayout1, blitPool1, blitSet[1]);

            // 你的代码：给第 0 帧建的卡槽
            VkDescriptorSetLayout denoiseLayout0, denoiseLayout1; VkDescriptorPool denoisePool0, denoisePool1; VkDescriptorSet denoiseSet[2];
            Lizeral::VulkanDescriptorBuilder()
                .BindImage(0, &denoiseInfos[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(1, &denoiseInfos[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(2, &denoiseInfos[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(3, &denoiseInfos[3], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(4, &denoiseInfos[4], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .Build(&vulkanDevice, denoiseLayout0, denoisePool0, denoiseSet[0]);

            Lizeral::VulkanDescriptorBuilder()
                .BindImage(0, &denoiseInfos[0], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(1, &denoiseInfos[1], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(2, &denoiseInfos[2], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(3, &denoiseInfos[3], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .BindImage(4, &denoiseInfos[4], VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT)
                .Build(&vulkanDevice, denoiseLayout1, denoisePool1, denoiseSet[1]);

            // =======================================================
            // ★ 构建所有 Pipeline
            // =======================================================
            VkShaderModule meshShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_mesh.spv"));
            VkShaderModule fragShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_frag.spv"));
            VkShaderModule taskShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_task.spv"));

            VkPushConstantRange pushRange{}; pushRange.stageFlags = VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT; pushRange.size = sizeof(PushConstants); 
            VkPipelineLayoutCreateInfo pipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; pipeLayoutInfo.pushConstantRangeCount = 1; pipeLayoutInfo.pPushConstantRanges = &pushRange; pipeLayoutInfo.setLayoutCount = 1; pipeLayoutInfo.pSetLayouts = &descriptorSetLayout;
            VkPipelineLayout pipelineLayout; vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &pipeLayoutInfo, nullptr, &pipelineLayout);

            VkPipeline graphicsPipeline = VulkanPipelineBuilder().AddShaderStage(VK_SHADER_STAGE_TASK_BIT_EXT, taskShaderModule).AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShaderModule).AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragShaderModule).SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST).SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, VK_FRONT_FACE_COUNTER_CLOCKWISE).SetMultisampling(VK_SAMPLE_COUNT_1_BIT).SetDepthStencil(true, true, VK_COMPARE_OP_GREATER_OR_EQUAL).DisableColorBlendAttachments(3).SetPipelineLayout(pipelineLayout).Build(&vulkanDevice, { VK_FORMAT_R8G8B8A8_UNORM, VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16_SFLOAT }, VK_FORMAT_D32_SFLOAT); 

            VkShaderModule lightVertModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "lighting_vert.spv"));
            VkShaderModule lightFragModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "lighting_frag.spv"));
            VkPushConstantRange lightPushRange{}; lightPushRange.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT; lightPushRange.size = sizeof(LightingPushConstants);
            VkDescriptorSetLayout lightSetLayouts[2] = { lightingLayouts[0], descriptorSetLayout };
            VkPipelineLayoutCreateInfo lightPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; lightPipeLayoutInfo.setLayoutCount = 2; lightPipeLayoutInfo.pSetLayouts = lightSetLayouts; lightPipeLayoutInfo.pushConstantRangeCount = 1; lightPipeLayoutInfo.pPushConstantRanges = &lightPushRange; 
            VkPipelineLayout lightingPipelineLayout; vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &lightPipeLayoutInfo, nullptr, &lightingPipelineLayout);
            VkPipeline lightingPipeline = VulkanPipelineBuilder().AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, lightVertModule).AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, lightFragModule).SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST).SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE).SetMultisampling(VK_SAMPLE_COUNT_1_BIT).SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS).AddColorBlendAttachment(false).AddColorBlendAttachment(false).SetPipelineLayout(lightingPipelineLayout).Build(&vulkanDevice, { VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED);

            VkShaderModule taaFragModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "taa_frag.spv"));
            VkPipelineLayoutCreateInfo taaPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; taaPipeLayoutInfo.setLayoutCount = 1; taaPipeLayoutInfo.pSetLayouts = &taaLayout0; 
            VkPipelineLayout taaPipelineLayout; vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &taaPipeLayoutInfo, nullptr, &taaPipelineLayout);
            VkPipeline taaPipeline = VulkanPipelineBuilder().AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, lightVertModule).AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, taaFragModule).SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST).SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE).SetMultisampling(VK_SAMPLE_COUNT_1_BIT).SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS).AddColorBlendAttachment(false).SetPipelineLayout(taaPipelineLayout).Build(&vulkanDevice, { VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED); 

            VkShaderModule blitFragModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "blit_frag.spv"));
            VkPipelineLayoutCreateInfo blitPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; blitPipeLayoutInfo.setLayoutCount = 1; blitPipeLayoutInfo.pSetLayouts = &blitLayout0; 
            VkPipelineLayout blitPipelineLayout; vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &blitPipeLayoutInfo, nullptr, &blitPipelineLayout);
            VkPipeline blitPipeline = VulkanPipelineBuilder().AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, lightVertModule).AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, blitFragModule).SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST).SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE).SetMultisampling(VK_SAMPLE_COUNT_1_BIT).SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS).AddColorBlendAttachment(false).SetPipelineLayout(blitPipelineLayout).Build(&vulkanDevice, { renderer.GetSwapchainFormat() }, VK_FORMAT_D32_SFLOAT);

            VkShaderModule denoiseFragModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "denoise_frag.spv"));
            VkPipelineLayoutCreateInfo denoisePipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
            denoisePipeLayoutInfo.setLayoutCount = 1; 
            denoisePipeLayoutInfo.pSetLayouts = &denoiseLayout0; 
            denoisePipeLayoutInfo.pushConstantRangeCount = 1; 
            denoisePipeLayoutInfo.pPushConstantRanges = &lightPushRange;

            VkPipelineLayout denoisePipelineLayout; 
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &denoisePipeLayoutInfo, nullptr, &denoisePipelineLayout);
            VkPipeline denoisePipeline = VulkanPipelineBuilder().AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, lightVertModule).AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, denoiseFragModule).SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST).SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE).SetMultisampling(VK_SAMPLE_COUNT_1_BIT).SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS).AddColorBlendAttachment(false).SetPipelineLayout(denoisePipelineLayout).Build(&vulkanDevice, { VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED); // 输出格式匹配 gDenoisedGI

            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), lightVertModule, nullptr); vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), lightFragModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), fragShaderModule, nullptr); vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), meshShaderModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), taskShaderModule, nullptr); vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), taaFragModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), blitFragModule, nullptr);  vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), denoiseFragModule, nullptr);

            auto CmdDrawMeshTasksEXT = (PFN_vkCmdDrawMeshTasksEXT)vkGetDeviceProcAddr(vulkanDevice.GetNativeDevice(), "vkCmdDrawMeshTasksEXT");

            bool firstFrame = true;
            Matrix4x4 prevViewProj; std::unordered_map<Entity, Matrix4x4> prevModelMats; bool isFirstFrameRun = true;
            static uint32_t frameIndex = 0; float lastTime = glfwGetTime();

            // 预分配 1000 个物体的台账空间
            std::vector<RTInstanceDesc> dummyInstances(1000); 
            VkBuffer rtInstBuffer;
            VkDeviceMemory rtInstMemory;
            uint64_t rtInstAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), dummyInstances, rtInstBuffer, rtInstMemory);
            g_AllocatedBuffers.push_back({rtInstBuffer, rtInstMemory}); // 自动回收

            std::cout << "\n[Sandbox] Dynamic Rendering & Ray Tracing Ready! Hold RMB + WASD to move." << std::endl;

            // ====================================================================
            // ★ 主循环
            // ====================================================================
            while (!glfwWindowShouldClose(window)) {
                
                float currentTime = glfwGetTime(); float dt = currentTime - lastTime; lastTime = currentTime;
                frameIndex++;

                uint32_t jitterIndex = frameIndex % 8 + 1; 
                float jitterX = CreateHaltonSequence(jitterIndex, 2) - 0.5f; float jitterY = CreateHaltonSequence(jitterIndex, 3) - 0.5f;
                float ndcJitterX = (jitterX * 2.0f) / static_cast<float>(WIDTH); float ndcJitterY = (jitterY * 2.0f) / static_cast<float>(HEIGHT);

                Input::GetInstance().Tick(); 
                if (Input::GetInstance().GetKey(Key::ESC)) glfwSetWindowShouldClose(window, true);
                glfwPollEvents();

                cameraControlSystem.Tick(dt, registry); cameraSystem.Tick(registry); 

                Matrix4x4 viewMat, projMatUnjittered, projMatJittered;
                auto cameraView = registry.view<TransformComponent, CameraComponent>();
                for (auto entity : cameraView) {
                    auto& camera = cameraView.get<CameraComponent>(entity);
                    camera.setProjectionMatrix(camera.BuildPerspectiveInfiniteReverseZ(camera.getFov(),camera.getAspect(),camera.getzNear()));
                    viewMat = camera.getViewMatrix(); 
                    projMatUnjittered = camera.BuildPerspectiveInfiniteReverseZ(camera.getFov(), camera.getAspect(), camera.getzNear());
                    projMatUnjittered[1][1] *= -1.0f; 
                    projMatJittered = projMatUnjittered;
                    projMatJittered[0][2] += ndcJitterX; projMatJittered[1][2] += ndcJitterY;
                    break; 
                }

                uint32_t ping = frameIndex % 2;       
                VkCommandBuffer cmd = renderer.BeginFrame();

                if (cmd != VK_NULL_HANDLE) {

                    // 1. 收集实例
                    RTInstanceDesc* mappedDesc;
                    vkMapMemory(vulkanDevice.GetNativeDevice(), rtInstMemory, 0, VK_WHOLE_SIZE, 0, (void**)&mappedDesc);

                    std::vector<VkAccelerationStructureInstanceKHR> tlasInstances;
                    uint32_t customInstanceId = 0;

                    for (auto entity : view) { // 借用外面的 view
                        auto& transform = view.get<TransformComponent>(entity);
                        auto& modelComp = view.get<VulkanModelComponent>(entity);
                        
                        // 安全检查：是否有模型、是否有 BLAS
                        if (!modelComp.IsLoaded() || !g_ModelCache[modelComp.getVulkanModelPath()].blas) continue;

                        // ★ 关键：在第一个循环里就把 res 取出来！
                        const auto& res = g_ModelCache[modelComp.getVulkanModelPath()];

                        // ★ 写入光追专属台账
                        mappedDesc[customInstanceId].vertexBuffer = res.vAddr;
                        mappedDesc[customInstanceId].indexBuffer  = res.globalIAddr; // 注意：这是上一轮加的专供光追的全局索引！
                        mappedDesc[customInstanceId].materialBuffer = res.matAddr;
                        mappedDesc[customInstanceId].textureOffset  = res.textureOffset;

                        // 获取矩阵并转换
                        Matrix4x4 modelMat = transform.getMatrix(); 
                        VkTransformMatrixKHR vkTransform{};
                        vkTransform.matrix[0][0] = modelMat[0][0]; vkTransform.matrix[0][1] = modelMat[0][1]; vkTransform.matrix[0][2] = modelMat[0][2]; vkTransform.matrix[0][3] = modelMat[0][3];
                        vkTransform.matrix[1][0] = modelMat[1][0]; vkTransform.matrix[1][1] = modelMat[1][1]; vkTransform.matrix[1][2] = modelMat[1][2]; vkTransform.matrix[1][3] = modelMat[1][3];
                        vkTransform.matrix[2][0] = modelMat[2][0]; vkTransform.matrix[2][1] = modelMat[2][1]; vkTransform.matrix[2][2] = modelMat[2][2]; vkTransform.matrix[2][3] = modelMat[2][3];

                        // 组装 TLAS 实例
                        VkAccelerationStructureInstanceKHR instance{};
                        instance.transform = vkTransform; 
                        
                        // ★ 把 customInstanceId 赋给 TLAS 节点，这样 Shader 里才能靠 ID 查到上面写的台账！
                        instance.instanceCustomIndex = customInstanceId++; 
                        
                        instance.mask = 0xFF; 
                        instance.flags = VK_GEOMETRY_INSTANCE_TRIANGLE_FACING_CULL_DISABLE_BIT_KHR; 
                        
                        // 顺便用取出的 res 替换掉之前那串长代码
                        instance.accelerationStructureReference = res.blas->GetDeviceAddress();
                        tlasInstances.push_back(instance);
                    }

                    // ★ 循环结束后，立刻解除映射！此时数据正式推入显存！
                    vkUnmapMemory(vulkanDevice.GetNativeDevice(), rtInstMemory);

                    // 2. 动态建树！
                    if (!tlasInstances.empty()) {
                        tlas.Build(cmd, ping, tlasInstances);
                        VkMemoryBarrier memoryBarrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
                        memoryBarrier.srcAccessMask = VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR;
                        memoryBarrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
                        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
                    }

                    // ====================================================================
                    // ★ 极速热更新 Descriptor：无论 Handle 是否改变，直接用最新地址覆盖！
                    // ====================================================================
                    VkAccelerationStructureKHR currentTlas = tlas.GetHandle(ping);
                    VkWriteDescriptorSetAccelerationStructureKHR asWrite{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET_ACCELERATION_STRUCTURE_KHR};
                    asWrite.accelerationStructureCount = 1;
                    asWrite.pAccelerationStructures = &currentTlas;

                    VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
                    write.dstSet = lightingSets[ping];
                    write.dstBinding = 4; 
                    write.dstArrayElement = 0;
                    write.descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
                    write.descriptorCount = 1;
                    write.pNext = &asWrite;
                    // GPU 会在绑定前自动应用这层更新
                    vkUpdateDescriptorSets(vulkanDevice.GetNativeDevice(), 1, &write, 0, nullptr);
                    
                    // --- G-Buffer Pass ---
                    VkImageLayout currentLayout = firstFrame ? VK_IMAGE_LAYOUT_UNDEFINED : VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
                    transitionImageLayout(cmd, gAlbedoMetallic.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gNormalRoughness.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gDepth.image, currentLayout, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);
                    transitionImageLayout(cmd, gVelocity.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    firstFrame = false;

                    VkRenderingAttachmentInfo colorAttachments[3] = {};
                    colorAttachments[0].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO; colorAttachments[0].imageView = gAlbedoMetallic.view; colorAttachments[0].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; colorAttachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; colorAttachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE; colorAttachments[0].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};
                    colorAttachments[1].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO; colorAttachments[1].imageView = gNormalRoughness.view; colorAttachments[1].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; colorAttachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; colorAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_STORE; colorAttachments[1].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};
                    colorAttachments[2].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO; colorAttachments[2].imageView = gVelocity.view; colorAttachments[2].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; colorAttachments[2].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; colorAttachments[2].storeOp = VK_ATTACHMENT_STORE_OP_STORE; colorAttachments[2].clearValue.color = {0.0f, 0.0f, 0.0f, 0.0f};
                    VkRenderingAttachmentInfo depthAttachment{VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO}; depthAttachment.imageView = gDepth.view; depthAttachment.imageLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL; depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE; depthAttachment.clearValue.depthStencil = {0.0f, 0};

                    VkRenderingInfo renderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO}; renderInfo.renderArea = {{0, 0}, {WIDTH, HEIGHT}}; renderInfo.layerCount = 1; renderInfo.colorAttachmentCount = 3; renderInfo.pColorAttachments = colorAttachments; renderInfo.pDepthAttachment = &depthAttachment;
                    vkCmdBeginRendering(cmd, &renderInfo);

                    VkViewport viewport{0, 0, (float)WIDTH, (float)HEIGHT, 0, 1}; vkCmdSetViewport(cmd, 0, 1, &viewport);
                    VkRect2D scissor{{0, 0}, {WIDTH, HEIGHT}}; vkCmdSetScissor(cmd, 0, 1, &scissor);
                    
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, graphicsPipeline);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout, 0, 1, &descriptorSet, 0, nullptr);

                    Matrix4x4 currentViewProjUnjittered = projMatUnjittered * viewMat;
                    if (isFirstFrameRun) prevViewProj = currentViewProjUnjittered;
                    
                    for (auto entity : view) {
                        auto& transform = view.get<TransformComponent>(entity);
                        auto& modelComp = view.get<VulkanModelComponent>(entity);
                        if (!modelComp.IsLoaded()) continue;
                        const auto& res = g_ModelCache[modelComp.getVulkanModelPath()];
                        Matrix4x4 currentModel = transform.getMatrix();
                        Matrix4x4 prevModel = isFirstFrameRun ? currentModel : prevModelMats[entity];

                        PushConstants pushData{};
                        pushData.mvp = (currentViewProjUnjittered * currentModel).transpose(); pushData.model = transform.getMatrix().transpose(); pushData.prevMvp = (prevViewProj * prevModel).transpose();
                        pushData.vertexBuffer = res.vAddr; pushData.meshletBuffer = res.mAddr; pushData.indexBuffer = res.iAddr; pushData.boundsBuffer = res.bAddr; pushData.materialBuffer = res.matAddr;
                        pushData.totalMeshlets = res.totalMeshlets; pushData.textureOffset = res.textureOffset; pushData.jitter = Vector2(ndcJitterX, ndcJitterY);

                        vkCmdPushConstants(cmd, pipelineLayout, VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(PushConstants), &pushData);
                        prevModelMats[entity] = currentModel;
                        CmdDrawMeshTasksEXT(cmd, (res.totalMeshlets + 63) / 64, 1, 1);
                    }
                    vkCmdEndRendering(cmd);

                    prevViewProj = currentViewProjUnjittered; isFirstFrameRun = false;

                    static bool isFirstTAAFrame = true;
                    if (isFirstTAAFrame) {
                        // transitionImageLayout(cmd, gSceneColor.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                        transitionImageLayout(cmd, gHistory[0].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                        transitionImageLayout(cmd, gHistory[1].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                        isFirstTAAFrame = false;
                    }
                    
                    transitionImageLayout(cmd, gAlbedoMetallic.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gNormalRoughness.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gDepth.image, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);
                    transitionImageLayout(cmd, gVelocity.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

                    // =======================================================
                    // --- 1. Lighting Pass (生成直接光与噪点GI) ---
                    // =======================================================
                    // 准备画布：每次画新帧，我们不在乎里面的旧数据 (UNDEFINED)，直接转为写入模式 (COLOR_ATTACHMENT)
                    transitionImageLayout(cmd, gDirectLight.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gNoisyGI.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    
                    VkRenderingAttachmentInfo lightAttachments[2] = {};
                    lightAttachments[0].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO; lightAttachments[0].imageView = gDirectLight.view; lightAttachments[0].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; lightAttachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; lightAttachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE; lightAttachments[0].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};
                    lightAttachments[1].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO; lightAttachments[1].imageView = gNoisyGI.view; lightAttachments[1].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; lightAttachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; lightAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_STORE; lightAttachments[1].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};

                    VkRenderingInfo lightRenderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO}; 
                    lightRenderInfo.renderArea = {{0, 0}, {WIDTH, HEIGHT}}; 
                    lightRenderInfo.layerCount = 1; 
                    lightRenderInfo.colorAttachmentCount = 2;
                    lightRenderInfo.pColorAttachments = lightAttachments;
                    
                    vkCmdBeginRendering(cmd, &lightRenderInfo); 
                    vkCmdSetViewport(cmd, 0, 1, &viewport); vkCmdSetScissor(cmd, 0, 1, &scissor);

                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, lightingPipeline);
                    VkDescriptorSet boundSets[2] = { lightingSets[ping], descriptorSet };
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, lightingPipelineLayout, 0, 2, boundSets, 0, nullptr);

                    LightingPushConstants lightPc{};
                    Matrix4x4 vpJittered = projMatJittered * viewMat;
                    lightPc.invViewProj = vpJittered.inverse().transpose();
                    lightPc.viewProj = vpJittered.transpose();
                    lightPc.cameraPos = registry.get<TransformComponent>(cameraEntity).getPosition(); 
                    lightPc.padding = 0.0f;
                    lightPc.frameIndex = frameIndex;
                    lightPc.padding2 = 0.0f;
                    lightPc.instanceDescAddr = rtInstAddr;
                    // lightPc.lightDir = Lizeral::Vector3(1.0f, 0.5f, 1.0f).normalize();
                    float time = glfwGetTime() * 0.2f; // 0.2 是旋转速度
                    lightPc.lightDir = Lizeral::Vector3(cos(time), 0.5f, sin(time)).normalize(); 
                    lightPc.lightPadding = 0.0f;
                    // 例如，给一个清晨/黄昏的暖橘色
                    lightPc.lightColor = Lizeral::Vector3(1.0f, 0.85f, 0.7f); 
                    
                    // 光强倍增器 (比如调到 4.0 甚至 10.0，让阳光刺眼)
                    lightPc.lightIntensity = 4.0f;


                    vkCmdPushConstants(cmd, lightingPipelineLayout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(LightingPushConstants), &lightPc);
                    vkCmdDraw(cmd, 3, 1, 0, 0);
                    vkCmdEndRendering(cmd);

                    // =======================================================
                    // --- 2. Denoise Pass (降噪并合成画面) ---
                    // =======================================================
                    // 1. 将刚画好的两张光照图转为读取模式 (READ_ONLY)，供给降噪器采样
                    transitionImageLayout(cmd, gDirectLight.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gNoisyGI.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    
                    // 2. 准备降噪合成后的输出画布
                    transitionImageLayout(cmd, gDenoisedGI.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

                    VkRenderingAttachmentInfo denoiseAttachment{VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO}; 
                    denoiseAttachment.imageView = gDenoisedGI.view; 
                    denoiseAttachment.imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; 
                    denoiseAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE; 
                    denoiseAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;

                    VkRenderingInfo denoiseRenderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO}; 
                    denoiseRenderInfo.renderArea = {{0, 0}, {WIDTH, HEIGHT}}; 
                    denoiseRenderInfo.layerCount = 1; 
                    denoiseRenderInfo.colorAttachmentCount = 1; 
                    denoiseRenderInfo.pColorAttachments = &denoiseAttachment;

                    vkCmdBeginRendering(cmd, &denoiseRenderInfo);
                    vkCmdSetViewport(cmd, 0, 1, &viewport); vkCmdSetScissor(cmd, 0, 1, &scissor);
                    
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, denoisePipeline);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, denoisePipelineLayout, 0, 1, &denoiseSet[ping], 0, nullptr);
                    vkCmdPushConstants(cmd, denoisePipelineLayout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(LightingPushConstants), &lightPc);
                    vkCmdDraw(cmd, 3, 1, 0, 0);
                    vkCmdEndRendering(cmd);

                    // =======================================================
                    // --- 3. TAA Pass (时空抗锯齿) ---
                    // =======================================================
                    // 1. 将刚合成好的画面转为读取模式，供给 TAA 采样
                    transitionImageLayout(cmd, gDenoisedGI.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    
                    // 2. 准备当前帧的 TAA 历史写入画布
                    transitionImageLayout(cmd, gHistory[ping].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    
                    VkRenderingAttachmentInfo taaAttachment{VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO}; 
                    taaAttachment.imageView = gHistory[ping].view; 
                    taaAttachment.imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; 
                    taaAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE; 
                    taaAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    
                    VkRenderingInfo taaRenderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO}; 
                    taaRenderInfo.renderArea = {{0, 0}, {WIDTH, HEIGHT}}; 
                    taaRenderInfo.layerCount = 1; 
                    taaRenderInfo.colorAttachmentCount = 1; 
                    taaRenderInfo.pColorAttachments = &taaAttachment;

                    vkCmdBeginRendering(cmd, &taaRenderInfo);
                    vkCmdSetViewport(cmd, 0, 1, &viewport); vkCmdSetScissor(cmd, 0, 1, &scissor);
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, taaPipeline);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, taaPipelineLayout, 0, 1, &taaSet[ping], 0, nullptr);
                    vkCmdDraw(cmd, 3, 1, 0, 0);
                    vkCmdEndRendering(cmd);

                    // =======================================================
                    // --- 4. Blit Pass (推上屏幕) ---
                    // =======================================================
                    // 将刚刚做完抗锯齿的画面，转为读取模式，送给屏幕
                    transitionImageLayout(cmd, gHistory[ping].image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    
                    renderer.BeginRendering(cmd); 
                    vkCmdSetViewport(cmd, 0, 1, &viewport); vkCmdSetScissor(cmd, 0, 1, &scissor);
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, blitPipeline);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, blitPipelineLayout, 0, 1, &blitSet[ping], 0, nullptr);
                    vkCmdDraw(cmd, 3, 1, 0, 0);

                    renderer.EndRendering(cmd);
                    renderer.EndFrame();
                }
            }

            std::cout << "--- EXITING PROGRAM ---" << std::endl;
            vkDeviceWaitIdle(vulkanDevice.GetNativeDevice());

            std::cout<<"[Pipeline] start destroying Pipeline"<<std::endl;
            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), blitPipeline, nullptr); vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), blitPipelineLayout, nullptr);
            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), taaPipeline, nullptr); vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), taaPipelineLayout, nullptr);
            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), lightingPipeline, nullptr); vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), lightingPipelineLayout, nullptr);
            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), graphicsPipeline, nullptr); vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), pipelineLayout, nullptr);
            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), denoisePipeline, nullptr); vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), denoisePipelineLayout, nullptr);

            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), descriptorPool, nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), lightPools[0], nullptr); vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), lightPools[1], nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), taaPool0, nullptr); vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), taaPool1, nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), blitPool0, nullptr); vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), blitPool1, nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), denoisePool0, nullptr); vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), denoisePool1, nullptr);

            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), descriptorSetLayout, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), lightingLayouts[0], nullptr); vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), lightingLayouts[1], nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), blitLayout0, nullptr); vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), blitLayout1, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), taaLayout0, nullptr); vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), taaLayout1, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), denoiseLayout0, nullptr); vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(), denoiseLayout1, nullptr); 

            vkDestroySampler(vulkanDevice.GetNativeDevice(), gBufferSampler, nullptr);

            std::cout<< "[VMA Debug] start clear BDA buffer"<<std::endl;
            for (auto& alloc : g_AllocatedBuffers) {
                vkDestroyBuffer(vulkanDevice.GetNativeDevice(), alloc.buffer, nullptr);
                vkFreeMemory(vulkanDevice.GetNativeDevice(), alloc.memory, nullptr);
            }

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gAlbedoMetallic.view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gAlbedoMetallic.image, gAlbedoMetallic.allocation);
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gNormalRoughness.view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gNormalRoughness.image, gNormalRoughness.allocation);
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gDepth.view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gDepth.image, gDepth.allocation);
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gVelocity.view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gVelocity.image, gVelocity.allocation);
            // vkDestroyImageView(vulkanDevice.GetNativeDevice(), gSceneColor.view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gSceneColor.image, gSceneColor.allocation);
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gHistory[0].view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gHistory[0].image, gHistory[0].allocation);
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gHistory[1].view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gHistory[1].image, gHistory[1].allocation);
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gDirectLight.view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gDirectLight.image, gDirectLight.allocation);
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gNoisyGI.view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gNoisyGI.image, gNoisyGI.allocation);
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gDenoisedGI.view, nullptr); vmaDestroyImage(vulkanDevice.GetAllocator(), gDenoisedGI.image, gDenoisedGI.allocation);
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