#include <iostream>
#include <exception>
#include <fstream>
#include <vector>
#include <cstring>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

// --- Vulkan & Engine RHI ---
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"
#include "runtime/function/render/rhi/vulkan/VulkanPipelineBuilder.h"
#include "runtime/function/render/MeshletBuilder/MeshletModelBuilder.h" 
#include "runtime/function/render/rhi/vulkan/VulkanCommandPool.h"
#include "runtime/function/render/rhi/vulkan/VulkanTexture.h"

// --- ECS & Camera System ---
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/input/input.h"
#include "runtime/function/render/CameraControlSystem/CameraControlSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h" 
#include "runtime/core/math/matrix4.h"

using namespace Lizeral;

const uint32_t WIDTH = 1280;
const uint32_t HEIGHT = 720;
const std::string SHADER_DIR = "C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/shader/";
const std::string MODEL_PATH = "C:/Lizeral Engine/LizeralEngine0.0.1/asset/maserati.glb";

struct PushConstants {
    Matrix4x4 mvp; 
    uint64_t vertexBuffer;
    uint64_t meshletBuffer;
    uint64_t indexBuffer;
    uint64_t boundsBuffer; 
    uint64_t materialBuffer;
    uint32_t totalMeshlets;
};

struct GBufferAttachment {
                VkImage image = VK_NULL_HANDLE;
                VmaAllocation allocation = VK_NULL_HANDLE;
                VkImageView view = VK_NULL_HANDLE;
                VkFormat format;
            };

// GLFW 输入回调代理
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

// Shader 读取助手
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

// BDA 缓冲分配助手
template <typename T>
uint64_t createBDABuffer(VkDevice device, VkPhysicalDevice physicalDevice, const std::vector<T>& data, VkBuffer& outBuffer, VkDeviceMemory& outMemory) {
    VkBufferCreateInfo bufferInfo{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bufferInfo.size = data.size() * sizeof(T);
    bufferInfo.usage = VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
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

// ====================================================================
// 2. 主程序入口
// ====================================================================
int main() {
    GLFWwindow* window = nullptr;
    VkSurfaceKHR surface = VK_NULL_HANDLE;
    VkInstance instance = VK_NULL_HANDLE;

    try {
        // --- 初始化 GLFW ---
        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API); 
        glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
        window = glfwCreateWindow(WIDTH, HEIGHT, "Lizeral Engine - Mesh Shader Rendering", nullptr, nullptr);

        glfwSetKeyCallback(window, glfwKeyCallback);
        glfwSetCursorPosCallback(window, glfwCursorPosCallback);
        glfwSetMouseButtonCallback(window, glfwMouseButtonCallback);

        // --- 初始化 Vulkan 核心架构 ---
        VulkanContext vulkanContext;
        uint32_t glfwExtensionCount = 0;
        const char** glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);
        std::vector<const char*> requiredExtensions(glfwExtensions, glfwExtensions + glfwExtensionCount);
        vulkanContext.Initialize("Lizeral Sandbox", requiredExtensions);
        instance = (VkInstance)vulkanContext.GetNativeInstance();

        if (glfwCreateWindowSurface(instance, window, nullptr, &surface) != VK_SUCCESS) 
            throw std::runtime_error("Failed to create window surface!");

        // 引入作用域以控制 Device 的生命周期
        {
            VulkanDevice vulkanDevice(&vulkanContext, surface);
            VulkanRenderer renderer(&vulkanContext, &vulkanDevice, window);

            Lizeral::VulkanCommandPool resourceCommandPool(&vulkanDevice);
            MeshletModelBuilder builder;

            builder.LoadAndSliceModel(MODEL_PATH);

            const auto& allTexturesData = builder.GetAllTextures();
            if (allTexturesData.empty()) throw std::runtime_error("No textures found in model!");

            std::vector<std::unique_ptr<Lizeral::VulkanTexture>> globalTextures;
            std::vector<VkDescriptorImageInfo> imageInfos;

            std::cout << "[Sandbox] Building Bindless Texture Array (" << allTexturesData.size() << " textures)..." << std::endl;
            
            for (const auto& texData : allTexturesData) {
                // 批量创建 VulkanTexture
                auto texture = std::make_unique<Lizeral::VulkanTexture>(
                    &vulkanDevice, 
                    &resourceCommandPool, 
                    texData.data(), 
                    static_cast<int>(texData.size())
                );
                
                // 收集描述符图片信息，准备塞进大数组
                VkDescriptorImageInfo info{};
                info.imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
                info.imageView = texture->GetImageView();
                info.sampler = texture->GetSampler();
                imageInfos.push_back(info);

                globalTextures.push_back(std::move(texture));
            }

            // ====================================================================
            // ★ 新阶段 2：创建 Bindless 描述符集布局
            // ====================================================================
            VkDescriptorSetLayoutBinding samplerLayoutBinding{};
            samplerLayoutBinding.binding = 0;
            // 告诉显卡：这是一个大小为 N 的数组！
            // samplerLayoutBinding.descriptorCount = static_cast<uint32_t>(globalTextures.size()); 
            samplerLayoutBinding.descriptorCount = 1024; 

            samplerLayoutBinding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            samplerLayoutBinding.pImmutableSamplers = nullptr;
            samplerLayoutBinding.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;

            // 【核心魔法】：声明这个数组是 Bindless 的，允许部分绑定和非均匀索引
            VkDescriptorBindingFlags bindlessFlags = VK_DESCRIPTOR_BINDING_PARTIALLY_BOUND_BIT;
            VkDescriptorSetLayoutBindingFlagsCreateInfo bindingFlagsInfo{};
            bindingFlagsInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_BINDING_FLAGS_CREATE_INFO;
            bindingFlagsInfo.bindingCount = 1;
            bindingFlagsInfo.pBindingFlags = &bindlessFlags;

            VkDescriptorSetLayoutCreateInfo layoutInfo{};
            layoutInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
            layoutInfo.bindingCount = 1;
            layoutInfo.pBindings = &samplerLayoutBinding;
            layoutInfo.pNext = &bindingFlagsInfo; // ★ 挂载 Bindless 标志！

            VkDescriptorSetLayout descriptorSetLayout;
            if (vkCreateDescriptorSetLayout(vulkanDevice.GetNativeDevice(), &layoutInfo, nullptr, &descriptorSetLayout) != VK_SUCCESS) {
                throw std::runtime_error("Failed to create bindless descriptor set layout!");
            }

            // ====================================================================
            // ★ 新阶段 3：创建描述符池与分配
            // ====================================================================
            VkDescriptorPoolSize poolSize{};
            poolSize.type = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            // poolSize.descriptorCount = static_cast<uint32_t>(globalTextures.size()); // 池子要足够大
            poolSize.descriptorCount = 1024; // 池子要足够大


            VkDescriptorPoolCreateInfo poolInfo{};
            poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
            poolInfo.poolSizeCount = 1;
            poolInfo.pPoolSizes = &poolSize;
            poolInfo.maxSets = 1;

            VkDescriptorPool descriptorPool;
            vkCreateDescriptorPool(vulkanDevice.GetNativeDevice(), &poolInfo, nullptr, &descriptorPool);

            VkDescriptorSetAllocateInfo allocInfo{};
            allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
            allocInfo.descriptorPool = descriptorPool;
            allocInfo.descriptorSetCount = 1;
            allocInfo.pSetLayouts = &descriptorSetLayout;

            VkDescriptorSet descriptorSet;
            vkAllocateDescriptorSets(vulkanDevice.GetNativeDevice(), &allocInfo, &descriptorSet);

            // ====================================================================
            // ★ 新阶段 4：将整个贴图数组，一次性砸进描述符集！
            // ====================================================================
            VkWriteDescriptorSet descriptorWrite{};
            descriptorWrite.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            descriptorWrite.dstSet = descriptorSet;
            descriptorWrite.dstBinding = 0;
            descriptorWrite.dstArrayElement = 0;
            descriptorWrite.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            // 告诉管线我要连续写多少张图进去
            descriptorWrite.descriptorCount = static_cast<uint32_t>(imageInfos.size()); 
            descriptorWrite.pImageInfo = imageInfos.data(); // 指向填满图片的数组数据指针

            vkUpdateDescriptorSets(vulkanDevice.GetNativeDevice(), 1, &descriptorWrite, 0, nullptr);
            
            auto CmdDrawMeshTasksEXT = (PFN_vkCmdDrawMeshTasksEXT)vkGetDeviceProcAddr(vulkanDevice.GetNativeDevice(), "vkCmdDrawMeshTasksEXT");

            VkBuffer vBuf, mBuf, iBuf, bBuf, matBuf;;
            VkDeviceMemory vMem, mMem, iMem, bMem, matMem;
            uint64_t vAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetVertices(), vBuf, vMem);
            uint64_t mAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetMeshlets(), mBuf, mMem);
            uint64_t iAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetMicroIndices(), iBuf, iMem);
            uint64_t bAddr = createBDABuffer(vulkanDevice.GetNativeDevice(), vulkanContext.GetPhysicalDevice(), builder.GetBounds(), bBuf, bMem);
            uint64_t matAddr = createBDABuffer(vulkanDevice.GetNativeDevice(),vulkanContext.GetPhysicalDevice(),builder.GetMaterials(),matBuf,matMem);

            auto createAttachment = [&](VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect) -> GBufferAttachment {
                GBufferAttachment attachment;
                attachment.format = format;

                // 1. 创建 Image (注意：WIDTH 和 HEIGHT 必须是你的窗口分辨率)
                VkImageCreateInfo imageInfo{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
                imageInfo.imageType = VK_IMAGE_TYPE_2D;
                imageInfo.extent.width = WIDTH;   // 你的全局窗口宽
                imageInfo.extent.height = HEIGHT; // 你的全局窗口高
                imageInfo.extent.depth = 1;
                imageInfo.mipLevels = 1;          // G-Buffer 绝对不需要 Mipmap！
                imageInfo.arrayLayers = 1;
                imageInfo.format = format;
                imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
                imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
                
                // ★ 延迟渲染的灵魂：既是渲染目标，又是可采样贴图！
                imageInfo.usage = usage | VK_IMAGE_USAGE_SAMPLED_BIT; 
                imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
                imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

                VmaAllocationCreateInfo allocInfo{};
                allocInfo.usage = VMA_MEMORY_USAGE_GPU_ONLY;

                if (vmaCreateImage(vulkanDevice.GetAllocator(), &imageInfo, &allocInfo, &attachment.image, &attachment.allocation, nullptr) != VK_SUCCESS) {
                    throw std::runtime_error("Failed to allocate G-Buffer image!");
                }

                // 2. 创建 Image View
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

            std::cout << "[Sandbox] Allocating G-Buffer Attachments..." << std::endl;

            // RT0: 漫反射颜色 (RGB) + 金属度 (A)
            GBufferAttachment gAlbedoMetallic = createAttachment(
                VK_FORMAT_R8G8B8A8_UNORM, 
                VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, 
                VK_IMAGE_ASPECT_COLOR_BIT
            );

            // RT1: 世界空间法线 (RGB) + 粗糙度 (A) -> 使用 16位浮点数保证法线精度！
            GBufferAttachment gNormalRoughness = createAttachment(
                VK_FORMAT_R16G16B16A16_SFLOAT, 
                VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, 
                VK_IMAGE_ASPECT_COLOR_BIT
            );

            // Depth: 深度图 (我们将用它在后续重构世界坐标，所以它也必须能被采样！)
            GBufferAttachment gDepth = createAttachment(
                VK_FORMAT_D32_SFLOAT, 
                VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT, 
                VK_IMAGE_ASPECT_DEPTH_BIT
            );

            VkSamplerCreateInfo gSamplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
            gSamplerInfo.magFilter = VK_FILTER_NEAREST;
            gSamplerInfo.minFilter = VK_FILTER_NEAREST;
            gSamplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
            gSamplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
            gSamplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
            gSamplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
            gSamplerInfo.maxAnisotropy = 1.0f;
            VkSampler gBufferSampler;
            vkCreateSampler(vulkanDevice.GetNativeDevice(), &gSamplerInfo, nullptr, &gBufferSampler);

            // 2. 布局：3 个 sampler2D
            std::vector<VkDescriptorSetLayoutBinding> lightingBindings(3);
            for(int i=0; i<3; i++) {
                lightingBindings[i].binding = i;
                lightingBindings[i].descriptorCount = 1;
                lightingBindings[i].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                lightingBindings[i].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }

            VkDescriptorSetLayoutCreateInfo lightLayoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            lightLayoutInfo.bindingCount = 3;
            lightLayoutInfo.pBindings = lightingBindings.data();
            VkDescriptorSetLayout lightingDescriptorSetLayout;
            vkCreateDescriptorSetLayout(vulkanDevice.GetNativeDevice(), &lightLayoutInfo, nullptr, &lightingDescriptorSetLayout);

            // 3. 分配描述符集 (★ 修复：为光照通道单独建一个小池子，因为之前的池子被 Bindless 用光了)
            VkDescriptorPoolSize lightPoolSize{};
            lightPoolSize.type = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            lightPoolSize.descriptorCount = 3; // 只需要装 3 张 G-Buffer 贴图

            VkDescriptorPoolCreateInfo lightPoolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            lightPoolInfo.poolSizeCount = 1;
            lightPoolInfo.pPoolSizes = &lightPoolSize;
            lightPoolInfo.maxSets = 1;

            VkDescriptorPool lightDescriptorPool;
            vkCreateDescriptorPool(vulkanDevice.GetNativeDevice(), &lightPoolInfo, nullptr, &lightDescriptorPool);

            VkDescriptorSetAllocateInfo lightAllocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            lightAllocInfo.descriptorPool = lightDescriptorPool; // ★ 使用全新的专属小池子
            lightAllocInfo.descriptorSetCount = 1;
            lightAllocInfo.pSetLayouts = &lightingDescriptorSetLayout;
            
            VkDescriptorSet lightingDescriptorSet;
            if (vkAllocateDescriptorSets(vulkanDevice.GetNativeDevice(), &lightAllocInfo, &lightingDescriptorSet) != VK_SUCCESS) {
                throw std::runtime_error("Failed to allocate lighting descriptor set!");
            }

            // 4. 写入 3 张 G-Buffer 贴图
            VkDescriptorImageInfo gInfos[3] = {};
            gInfos[0].imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
            gInfos[0].imageView = gAlbedoMetallic.view; gInfos[0].sampler = gBufferSampler;
            
            gInfos[1].imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
            gInfos[1].imageView = gNormalRoughness.view; gInfos[1].sampler = gBufferSampler;
            
            gInfos[2].imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
            gInfos[2].imageView = gDepth.view; gInfos[2].sampler = gBufferSampler;

            std::vector<VkWriteDescriptorSet> lightWrites(3);
            for(int i=0; i<3; i++) {
                lightWrites[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                lightWrites[i].dstSet = lightingDescriptorSet;
                lightWrites[i].dstBinding = i;
                lightWrites[i].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                lightWrites[i].descriptorCount = 1;
                lightWrites[i].pImageInfo = &gInfos[i];
            }
            vkUpdateDescriptorSets(vulkanDevice.GetNativeDevice(), 3, lightWrites.data(), 0, nullptr);

            // --- 管线构建 ---
            VkShaderModule meshShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_mesh.spv"));
            VkShaderModule fragShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_frag.spv"));
            VkShaderModule taskShaderModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "car_task.spv"));

            VkPushConstantRange pushConstantRange{};
            pushConstantRange.stageFlags = VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT; // 当前仅用 Mesh
            pushConstantRange.offset = 0;
            pushConstantRange.size = sizeof(PushConstants); 

            VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
            pipelineLayoutInfo.pushConstantRangeCount = 1;
            pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;

            pipelineLayoutInfo.setLayoutCount = 1; 
            pipelineLayoutInfo.pSetLayouts = &descriptorSetLayout;
            
            VkPipelineLayout pipelineLayout;
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &pipelineLayoutInfo, nullptr, &pipelineLayout);

            VulkanPipelineBuilder builderPipeline;
            VkPipeline graphicsPipeline = builderPipeline
                .AddShaderStage(VK_SHADER_STAGE_TASK_BIT_EXT, taskShaderModule)
                .AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShaderModule) 
                .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragShaderModule)
                // 明确声明：我要画三角形
                .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
                // 明确声明：填充模式，不剔除背面
                .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, VK_FRONT_FACE_COUNTER_CLOCKWISE)
                // 明确声明：不抗锯齿
                .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
                // 明确声明：开启深度测试，允许写入深度
                .SetDepthStencil(true, true, VK_COMPARE_OP_LESS)
                // ★ 明确声明：我的 G-Buffer 有 2 个 Render Target，它们都绝对不需要颜色混合！
                .DisableColorBlendAttachments(2)
                .SetPipelineLayout(pipelineLayout)
                // ★ 传入 MRT 格式数组！(注意：这里的格式必须和我们将要创建的 VkImage 格式一模一样)
                .Build(&vulkanDevice, 
                       { VK_FORMAT_R8G8B8A8_UNORM, VK_FORMAT_R16G16B16A16_SFLOAT }, // RT0, RT1
                       VK_FORMAT_D32_SFLOAT); // Depth

            // ====================================================================
            // ★ 新阶段：构建 Lighting Pipeline
            // ====================================================================
            VkShaderModule lightVertModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "lighting_vert.spv"));
            VkShaderModule lightFragModule = createShaderModule(vulkanDevice.GetNativeDevice(), readFile(SHADER_DIR + "lighting_frag.spv"));

            VkPipelineLayoutCreateInfo lightPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
            lightPipeLayoutInfo.setLayoutCount = 1;
            lightPipeLayoutInfo.pSetLayouts = &lightingDescriptorSetLayout; // 绑定刚才建的 Layout
            VkPipelineLayout lightingPipelineLayout;
            vkCreatePipelineLayout(vulkanDevice.GetNativeDevice(), &lightPipeLayoutInfo, nullptr, &lightingPipelineLayout);

            VulkanPipelineBuilder builderLighting;
            VkPipeline lightingPipeline = builderLighting
                .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, lightVertModule)
                .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, lightFragModule)
                .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
                .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE) // 剔除其实无所谓，因为是全屏 2D
                .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
                .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS) // ★ 光照阶段绝对不需要深度测试！
                .AddColorBlendAttachment(false) // 直接覆盖 Swapchain 像素
                .SetPipelineLayout(lightingPipelineLayout)
                .Build(&vulkanDevice, { renderer.GetSwapchainFormat() }, VK_FORMAT_D32_SFLOAT); // ★ 目标是屏幕，不需要深度图！


            //destroyShaderModule
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), lightVertModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), lightFragModule, nullptr);

            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), fragShaderModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), meshShaderModule, nullptr);
            vkDestroyShaderModule(vulkanDevice.GetNativeDevice(), taskShaderModule, nullptr);

            // --- 场景 ECS 装配 ---
            Registry registry;
            CameraControlSystem cameraControlSystem;
            CameraSystem cameraSystem;

            Entity cameraEntity = registry.create();
            {
                auto& t = registry.emplace<TransformComponent>(cameraEntity);
                t.setPosition(Vector3(0.0f, 0.0f, 1.0f)); 
                
                auto& cam = registry.emplace<CameraComponent>(cameraEntity);
                cam.setFov(45.0f); cam.setAspect((float)WIDTH / (float)HEIGHT);
                cam.setzNear(0.001f); cam.setzFar(1000.0f); cam.setMain(true);

                auto& ctrl = registry.emplace<CameraControlComponent>(cameraEntity);
                ctrl.setMoveSpeed(1.0f); ctrl.setSensitivityX(0.1f); ctrl.setSensitivityY(0.1f);
                ctrl.setYaw(0.0f);
            }

            std::cout << "\n[Sandbox] Dynamic Rendering & ECS Ready! Hold RMB + WASD to move." << std::endl;

            auto transitionImageLayout = [&](VkCommandBuffer cmd, VkImage image, VkImageLayout oldLayout, VkImageLayout newLayout, VkImageAspectFlags aspectMask) {
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
                } else {
                    throw std::invalid_argument("Unsupported layout transition!");
                }

                vkCmdPipelineBarrier(cmd, sourceStage, destinationStage, 0, 0, nullptr, 0, nullptr, 1, &barrier);
            };
            
            // --- 渲染主循环 ---
            float lastTime = glfwGetTime();
            while (!glfwWindowShouldClose(window)) {
                
                float currentTime = glfwGetTime();
                float dt = currentTime - lastTime;
                lastTime = currentTime;

                Input::GetInstance().Tick(); 
                if (Input::GetInstance().GetKey(Key::ESC)) glfwSetWindowShouldClose(window, true);
                glfwPollEvents();

                cameraControlSystem.Tick(dt, registry);
                cameraSystem.Tick(registry); 

                // 提取相机矩阵
                Matrix4x4 viewMat, projMat;
                auto view = registry.view<TransformComponent, CameraComponent>();
                for (auto entity : view) {
                    auto& camera = view.get<CameraComponent>(entity);
                    camera.setProjectionMatrix(camera.buildPerspective(camera.getFov(), camera.getAspect(), camera.getzNear(), camera.getzFar()));
                    viewMat = camera.getViewMatrix(); 
                    projMat = camera.getProjectionMatrix();
                    break; 
                }

                projMat[1][1] *= -1.0f; // Vulkan Y 轴补丁
                Matrix4x4 modelMat = Matrix4x4::getScale(Vector3(10.0f, 10.0f, 10.0f));

                PushConstants pushData{};
                pushData.mvp = (projMat * viewMat * modelMat).transpose(); 
                pushData.vertexBuffer = vAddr;
                pushData.meshletBuffer = mAddr;
                pushData.indexBuffer = iAddr;
                pushData.boundsBuffer = bAddr; 
                pushData.materialBuffer = matAddr;
                pushData.totalMeshlets = builder.GetMeshlets().size(); 

                // 提交渲染
                VkCommandBuffer cmd = renderer.BeginFrame();
                if (cmd != VK_NULL_HANDLE) {
                    // ====================================================================
                    // ★ 劫持阶段 1：布局转换 (准备写入 G-Buffer)
                    // ====================================================================
                    transitionImageLayout(cmd, gAlbedoMetallic.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gNormalRoughness.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gDepth.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);

                    // ====================================================================
                    // ★ 劫持阶段 2：配置自定义的 G-Buffer 渲染目标
                    // ====================================================================
                    VkRenderingAttachmentInfo colorAttachments[2] = {};
                    
                    // RT0: 清屏颜色设为纯黑
                    colorAttachments[0].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    colorAttachments[0].imageView = gAlbedoMetallic.view;
                    colorAttachments[0].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                    colorAttachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
                    colorAttachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    colorAttachments[0].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};

                    // RT1: 清屏法线为0，粗糙度默认为1.0
                    colorAttachments[1].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    colorAttachments[1].imageView = gNormalRoughness.view;
                    colorAttachments[1].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                    colorAttachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
                    colorAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    colorAttachments[1].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};

                    VkRenderingAttachmentInfo depthAttachment{};
                    depthAttachment.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                    depthAttachment.imageView = gDepth.view;
                    depthAttachment.imageLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
                    depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
                    depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                    depthAttachment.clearValue.depthStencil = {1.0f, 0};

                    VkRenderingInfo renderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO};
                    renderInfo.renderArea = {{0, 0}, {WIDTH, HEIGHT}};
                    renderInfo.layerCount = 1;
                    renderInfo.colorAttachmentCount = 2;
                    renderInfo.pColorAttachments = colorAttachments;
                    renderInfo.pDepthAttachment = &depthAttachment;

                    // ====================================================================
                    // ★ 劫持阶段 3：开始疯狂输出几何数据到 G-Buffer！
                    // ====================================================================
                    auto vkCmdBeginRendering = (PFN_vkCmdBeginRendering)vkGetDeviceProcAddr(vulkanDevice.GetNativeDevice(), "vkCmdBeginRendering");
                    auto vkCmdEndRendering = (PFN_vkCmdEndRendering)vkGetDeviceProcAddr(vulkanDevice.GetNativeDevice(), "vkCmdEndRendering");

                    vkCmdBeginRendering(cmd, &renderInfo);

                    VkViewport viewport{};
                    viewport.x = 0.0f;
                    viewport.y = 0.0f;
                    viewport.width = static_cast<float>(WIDTH);
                    viewport.height = static_cast<float>(HEIGHT);
                    viewport.minDepth = 0.0f;
                    viewport.maxDepth = 1.0f;
                    vkCmdSetViewport(cmd, 0, 1, &viewport);

                    VkRect2D scissor{};
                    scissor.offset = {0, 0};
                    scissor.extent = {static_cast<uint32_t>(WIDTH), static_cast<uint32_t>(HEIGHT)};
                    vkCmdSetScissor(cmd, 0, 1, &scissor);
                    
                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, graphicsPipeline);
                    vkCmdPushConstants(cmd, pipelineLayout, VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(PushConstants), &pushData);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout, 0, 1, &descriptorSet, 0, nullptr);
                    
                    uint32_t taskGroupCount = (builder.GetMeshlets().size() + 63) / 64;
                    CmdDrawMeshTasksEXT(cmd, taskGroupCount, 1, 1);

                    vkCmdEndRendering(cmd);

                    // ====================================================================
                    // ★ 劫持阶段 4：布局转换 (准备作为贴图给光照阶段采样！)
                    // ====================================================================
                    transitionImageLayout(cmd, gAlbedoMetallic.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gNormalRoughness.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                    transitionImageLayout(cmd, gDepth.image, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);

                    // --- 渲染真正的屏幕 (Lighting Pass) 将在这里进行 ---
                    // 现在屏幕什么都没画，但我们必须保证 Swapchain 的状态机正常闭环
                    renderer.BeginRendering(cmd); 
                    
                    // 重新设置视口和裁剪 (因为我们换了一个 Render Pass)
                    vkCmdSetViewport(cmd, 0, 1, &viewport);
                    vkCmdSetScissor(cmd, 0, 1, &scissor);

                    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, lightingPipeline);
                    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, lightingPipelineLayout, 0, 1, &lightingDescriptorSet, 0, nullptr);
                    
                    // ★ 神奇的发车指令：画 3 个顶点，不需要顶点 Buffer！
                    vkCmdDraw(cmd, 3, 1, 0, 0);

                    renderer.EndRendering(cmd);

                    renderer.EndFrame();
                }
            }
            
            vkDeviceWaitIdle(vulkanDevice.GetNativeDevice());

            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), lightingPipeline, nullptr);
            vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), lightingPipelineLayout, nullptr);

            vkDestroyPipeline(vulkanDevice.GetNativeDevice(), graphicsPipeline, nullptr);
            vkDestroyPipelineLayout(vulkanDevice.GetNativeDevice(), pipelineLayout, nullptr);

            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(),descriptorPool,nullptr);
            vkDestroyDescriptorPool(vulkanDevice.GetNativeDevice(), lightDescriptorPool, nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(),descriptorSetLayout,nullptr);

            vkDestroySampler(vulkanDevice.GetNativeDevice(),gBufferSampler,nullptr);
            vkDestroyDescriptorSetLayout(vulkanDevice.GetNativeDevice(),lightingDescriptorSetLayout,nullptr);
            
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), vBuf, nullptr); vkFreeMemory(vulkanDevice.GetNativeDevice(), vMem, nullptr);
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), mBuf, nullptr); vkFreeMemory(vulkanDevice.GetNativeDevice(), mMem, nullptr);
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), iBuf, nullptr); vkFreeMemory(vulkanDevice.GetNativeDevice(), iMem, nullptr);
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), bBuf, nullptr); vkFreeMemory(vulkanDevice.GetNativeDevice(), bMem, nullptr);
            vkDestroyBuffer(vulkanDevice.GetNativeDevice(), matBuf, nullptr); vkFreeMemory(vulkanDevice.GetNativeDevice(), matMem, nullptr);

            // 清理 G-Buffer
            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gAlbedoMetallic.view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gAlbedoMetallic.image, gAlbedoMetallic.allocation);

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gNormalRoughness.view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gNormalRoughness.image, gNormalRoughness.allocation);

            vkDestroyImageView(vulkanDevice.GetNativeDevice(), gDepth.view, nullptr);
            vmaDestroyImage(vulkanDevice.GetAllocator(), gDepth.image, gDepth.allocation);

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