#include "VulkanRenderingSystem.h"
#include "runtime/function/render/rhi/vulkan/VulkanPipelineBuilder.h"
#include "runtime/function/render/rhi/vulkan/VulkanDescriptorBuilder.h"
#include "runtime/function/render/MeshletBuilder/MeshletModelBuilder.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Model/VulkanModelComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Light/DirectionalLightComponent.h"

#include <iostream>
#include <fstream>

namespace Lizeral {

    // =========================================================================================
    // 辅助工具：利用你新写的 VulkanBuffer 一键创建 BDA 缓冲并写入数据
    // =========================================================================================

    template <typename T>
    std::unique_ptr<VulkanBuffer> VulkanRenderingSystem::CreateBDABuffer(const std::vector<T>& data) {
        if (data.empty()) return nullptr;

        VkDeviceSize bufferSize = data.size() * sizeof(T);
        
        // 使用你的 VulkanBuffer，指定 VMA_MEMORY_USAGE_CPU_TO_GPU 以便 WriteData
        auto buffer = std::make_unique<VulkanBuffer>(
            m_device, 
            bufferSize, 
            VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT | VK_BUFFER_USAGE_ACCELERATION_STRUCTURE_BUILD_INPUT_READ_ONLY_BIT_KHR,
            VMA_MEMORY_USAGE_CPU_TO_GPU 
        );
        
        buffer->WriteData(data.data(), bufferSize);
        return buffer;
    }

    GBufferAttachment VulkanRenderingSystem::CreateAttachment(VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect) {
        GBufferAttachment attachment; 
        attachment.format = format;

        VkImageCreateInfo imageInfo{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
        imageInfo.imageType = VK_IMAGE_TYPE_2D; 
        imageInfo.extent = {m_width, m_height, 1};
        imageInfo.mipLevels = 1; 
        imageInfo.arrayLayers = 1; 
        imageInfo.format = format;
        imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL; 
        imageInfo.usage = usage | VK_IMAGE_USAGE_SAMPLED_BIT; 
        imageInfo.samples = VK_SAMPLE_COUNT_1_BIT; 
        imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

        VmaAllocationCreateInfo allocInfo{}; 
        allocInfo.usage = VMA_MEMORY_USAGE_GPU_ONLY;
        
        vmaCreateImage(m_device->GetAllocator(), &imageInfo, &allocInfo, &attachment.image, &attachment.allocation, nullptr);

        VkImageViewCreateInfo viewInfo{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
        viewInfo.image = attachment.image; 
        viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D; 
        viewInfo.format = format;
        viewInfo.subresourceRange.aspectMask = aspect; 
        viewInfo.subresourceRange.levelCount = 1; 
        viewInfo.subresourceRange.layerCount = 1;
        
        vkCreateImageView(m_device->GetNativeDevice(), &viewInfo, nullptr, &attachment.view);
        return attachment;
    }

    void VulkanRenderingSystem::DestroyAttachment(GBufferAttachment& attachment) {
        if (attachment.view != VK_NULL_HANDLE) {
            vkDestroyImageView(m_device->GetNativeDevice(), attachment.view, nullptr);
            attachment.view = VK_NULL_HANDLE;
        }
        if (attachment.image != VK_NULL_HANDLE && attachment.allocation != VK_NULL_HANDLE) {
            vmaDestroyImage(m_device->GetAllocator(), attachment.image, attachment.allocation);
            attachment.image = VK_NULL_HANDLE;
            attachment.allocation = VK_NULL_HANDLE;
        }
    }

    float VulkanRenderingSystem::CreateHaltonSequence(uint32_t index, uint32_t base) {
        float f = 1.0f; float result = 0.0f; uint32_t i = index;
        while (i > 0) {
            f = f / static_cast<float>(base);
            result = result + f * static_cast<float>(i % base);
            i = static_cast<uint32_t>(std::floor(static_cast<float>(i) / static_cast<float>(base)));
        }
        return result;
    }

    void VulkanRenderingSystem::TransitionImageLayout(VkCommandBuffer cmd, VkImage image, VkImageLayout oldLayout, VkImageLayout newLayout, VkImageAspectFlags aspectMask) {
        VkImageMemoryBarrier barrier{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
        barrier.oldLayout = oldLayout; barrier.newLayout = newLayout;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED; barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = image; barrier.subresourceRange.aspectMask = aspectMask;
        barrier.subresourceRange.baseMipLevel = 0; barrier.subresourceRange.levelCount = 1;
        barrier.subresourceRange.baseArrayLayer = 0; barrier.subresourceRange.layerCount = 1;

        VkPipelineStageFlags sourceStage; VkPipelineStageFlags destinationStage;

        if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL) {
            barrier.srcAccessMask = 0; barrier.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
            sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT; destinationStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        } else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
            barrier.srcAccessMask = 0; barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
            sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT; destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
        } else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL) {
            barrier.srcAccessMask = 0; barrier.dstAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
            sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT; destinationStage = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT;
        } else if (oldLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
            barrier.srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT; barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
            sourceStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT; destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
        } else if (oldLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
            barrier.srcAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT; barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
            sourceStage = VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT; destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
        } else if (oldLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL) {
            barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT; barrier.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
            sourceStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT; destinationStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        } 
        else if (oldLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL) {
            barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT;
            barrier.dstAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
            sourceStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
            destinationStage = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT; 
        } 
        else {
            std::cerr << "[VulkanRenderingSystem] FATAL: Unsupported layout transition from " << oldLayout << " to " << newLayout << std::endl;
            throw std::invalid_argument("Unsupported layout transition!");
        }

        vkCmdPipelineBarrier(cmd, sourceStage, destinationStage, 0, 0, nullptr, 0, nullptr, 1, &barrier);
    }


    void VulkanRenderingSystem::Initialize(VulkanContext* context, VulkanDevice* device, VulkanRenderer* renderer, uint32_t width, uint32_t height) {
        m_context = context;
        m_device = device;
        m_renderer = renderer;
        m_width = width;
        m_height = height;

        std::cout << "[VulkanRenderingSystem] Initializing..." << std::endl;

        struct RTInstanceDesc {
            uint64_t vertexBuffer; uint64_t indexBuffer; uint64_t materialBuffer;
            uint32_t textureOffset; uint32_t padding[3];
        };

        // 1. 创建基础池与光追结构
        if (!m_tlas) {
            m_resourceCommandPool = std::make_unique<VulkanCommandPool>(m_device);
            m_tlas = std::make_unique<VulkanTLAS>(m_device, 2);

            VulkanCommandBuffer tempCmd(m_device, m_resourceCommandPool.get());
            tempCmd.Begin(VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT);
            std::vector<VkAccelerationStructureInstanceKHR> emptyInstances;
            m_tlas->Build(tempCmd.GetNativeBuffer(), 0, emptyInstances);
            m_tlas->Build(tempCmd.GetNativeBuffer(), 1, emptyInstances);
            tempCmd.End();
            tempCmd.SubmitAndIdle();

            std::vector<RTInstanceDesc> dummyInstances(1000); 
            m_rtInstanceBuffer = CreateBDABuffer(dummyInstances);

            std::vector<GPUInstanceData> dummyGPU(10000);   
            m_globalInstanceBuffer = CreateBDABuffer(dummyGPU);
        }

        CreateAttachments();

        // 4. 创建全局采样器
        VkSamplerCreateInfo samplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
        samplerInfo.magFilter = VK_FILTER_LINEAR; 
        samplerInfo.minFilter = VK_FILTER_LINEAR;
        samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
        samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE; 
        samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE; 
        samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
        vkCreateSampler(m_device->GetNativeDevice(), &samplerInfo, nullptr, &m_gBufferSampler);

        // 5. 编译与构建管线
        BuildPipelines();

        m_firstFrame = true;
        m_isFirstFrameRun = true;

        std::cout << "[VulkanRenderingSystem] Initialization Complete." << std::endl;
    }

    void VulkanRenderingSystem::CreateAttachments() {

        m_gAlbedoMetallic = CreateAttachment(VK_FORMAT_R8G8B8A8_UNORM, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        m_gNormalRoughness = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        m_gDepth = CreateAttachment(VK_FORMAT_D32_SFLOAT, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT, VK_IMAGE_ASPECT_DEPTH_BIT);
        m_gVelocity = CreateAttachment(VK_FORMAT_R16G16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        
        m_gDirectLight = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        m_gNoisyGI = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        m_gDenoisedGI = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        m_gDenoisedGITemp = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT, VK_IMAGE_ASPECT_COLOR_BIT);

        m_gGIHistory[0] = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        m_gGIHistory[1] = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        
        m_gMomentsHistory[0] = CreateAttachment(VK_FORMAT_R16G16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        m_gMomentsHistory[1] = CreateAttachment(VK_FORMAT_R16G16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        
        m_gHistory[0] = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
        m_gHistory[1] = CreateAttachment(VK_FORMAT_R16G16B16A16_SFLOAT, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
    }


    void VulkanRenderingSystem::Shutdown() {
        if (!m_device) return;
        
        // 1. 第一步永远是死等 GPU 停下手中的所有活！
        vkDeviceWaitIdle(m_device->GetNativeDevice());
        VkDevice device = m_device->GetNativeDevice();

        std::cout << "[RenderingSystem] Shutting down..." << std::endl;

        // =======================================================
        // ★ 修复区 1：清理所有包含显存的容器
        // =======================================================
        m_modelCache.clear();      // 释放所有物体的 Vertex/Index BDA Buffer
        m_globalTextures.clear();  // 释放所有 Bindless 贴图

        // =======================================================
        // ★ 修复区 2：手动销毁所有动态创建的全局大 Buffer (核心漏网之鱼)
        // =======================================================
        if (m_rtInstanceBuffer) m_rtInstanceBuffer.reset();
        if (m_globalInstanceBuffer) m_globalInstanceBuffer.reset(); // 释放 10000 个矩阵的显存
        if (m_frameResource.buffer) m_frameResource.buffer.reset(); // 释放 FrameData 显存

        // =======================================================
        // ★ 修复区 3：销毁常规 Vulkan 资源
        // =======================================================
        DestroyPipelines(device);
        DestroyDescriptors(device);
        DestroyAttachments();

        if (m_gBufferSampler) {
            vkDestroySampler(device, m_gBufferSampler, nullptr);
            m_gBufferSampler = VK_NULL_HANDLE;
        }

        // =======================================================
        // ★ 修复区 4：销毁光追与命令池
        // =======================================================
        if (m_tlas) m_tlas.reset();
        if (m_resourceCommandPool) m_resourceCommandPool.reset();

        std::cout << "[RenderingSystem] Shutdown complete." << std::endl;
    }

    void VulkanRenderingSystem::DestroyPipelines(VkDevice device){
        if (m_blitPipeline)       vkDestroyPipeline(device, m_blitPipeline, nullptr);
        if (m_blitPipelineLayout) vkDestroyPipelineLayout(device, m_blitPipelineLayout, nullptr);

        if (m_taaPipeline)        vkDestroyPipeline(device, m_taaPipeline, nullptr);
        if (m_taaPipelineLayout)  vkDestroyPipelineLayout(device, m_taaPipelineLayout, nullptr);

        if (m_svgfTemporalPipeline)       vkDestroyPipeline(device, m_svgfTemporalPipeline, nullptr);
        if (m_svgfTemporalPipelineLayout) vkDestroyPipelineLayout(device, m_svgfTemporalPipelineLayout, nullptr);

        if (m_svgfATrousPipeline)       vkDestroyPipeline(device, m_svgfATrousPipeline, nullptr);
        if (m_svgfATrousPipelineLayout) vkDestroyPipelineLayout(device, m_svgfATrousPipelineLayout, nullptr);

        if (m_lightingPipeline)       vkDestroyPipeline(device, m_lightingPipeline, nullptr);
        if (m_lightingPipelineLayout) vkDestroyPipelineLayout(device, m_lightingPipelineLayout, nullptr);

        if (m_graphicsPipeline)       vkDestroyPipeline(device, m_graphicsPipeline, nullptr);
        if (m_graphicsPipelineLayout) vkDestroyPipelineLayout(device, m_graphicsPipelineLayout, nullptr);
    }

    void VulkanRenderingSystem::DestroyDescriptors(VkDevice device){
        if (m_descriptorPool) vkDestroyDescriptorPool(device, m_descriptorPool, nullptr);
        for (int i = 0; i < 2; i++) {
            if (m_lightPools[i])   vkDestroyDescriptorPool(device, m_lightPools[i], nullptr);
            if (m_svgfTemporalPools[i]) vkDestroyDescriptorPool(device, m_svgfTemporalPools[i], nullptr);
            if (m_svgfATrousPools[i]) vkDestroyDescriptorPool(device, m_svgfATrousPools[i], nullptr);
            if (m_taaPools[i])     vkDestroyDescriptorPool(device, m_taaPools[i], nullptr);
            if (m_blitPools[i])    vkDestroyDescriptorPool(device, m_blitPools[i], nullptr);
        }

        if (m_descriptorSetLayout) vkDestroyDescriptorSetLayout(device, m_descriptorSetLayout, nullptr);
        for (int i = 0; i < 2; i++) {
            if (m_lightingLayouts[i]) vkDestroyDescriptorSetLayout(device, m_lightingLayouts[i], nullptr);
            if (m_svgfTemporalLayouts[i])  vkDestroyDescriptorSetLayout(device, m_svgfTemporalLayouts[i], nullptr);
            if (m_svgfATrousLayouts[i])  vkDestroyDescriptorSetLayout(device, m_svgfATrousLayouts[i], nullptr);
            if (m_taaLayouts[i])      vkDestroyDescriptorSetLayout(device, m_taaLayouts[i], nullptr);
            if (m_blitLayouts[i])     vkDestroyDescriptorSetLayout(device, m_blitLayouts[i], nullptr);
        }
    }

    void VulkanRenderingSystem::DestroyAttachments(){
        DestroyAttachment(m_gAlbedoMetallic); 
        DestroyAttachment(m_gNormalRoughness);
        DestroyAttachment(m_gDepth);          
        DestroyAttachment(m_gVelocity);
        DestroyAttachment(m_gDirectLight);    
        DestroyAttachment(m_gNoisyGI);
        DestroyAttachment(m_gDenoisedGI);  
        DestroyAttachment(m_gDenoisedGITemp);
        DestroyAttachment(m_gGIHistory[0]);
        DestroyAttachment(m_gGIHistory[1]); 
        DestroyAttachment(m_gMomentsHistory[0]);
        DestroyAttachment(m_gMomentsHistory[1]);
        DestroyAttachment(m_gHistory[0]); 
        DestroyAttachment(m_gHistory[1]);
    }


    void VulkanRenderingSystem::Resize(uint32_t width, uint32_t height) {
        if (width == 0 || height == 0) return;
        vkDeviceWaitIdle(m_device->GetNativeDevice());
        
        // 只更新画板（G-Buffer）真实尺寸
        m_width = width;
        m_height = height;

        DestroyAttachments();
        CreateAttachments();

        UpdateDescriptorSets();

        m_isFirstFrameRun = true;
        m_firstFrame = true;
    }

    void VulkanRenderingSystem::SetViewport(int x, int y, uint32_t width, uint32_t height) {
        m_viewX = x;
        m_viewY = y;
        m_viewW = width;
        m_viewH = height;
    }

    VulkanModelResource& VulkanRenderingSystem::GetOrLoadModel(const std::string& path) {
        // 如果已经加载过，直接返回缓存中的引用
        if (m_modelCache.find(path) != m_modelCache.end()) {
            return m_modelCache[path];
        }

        std::cout << "[RenderingSystem] Loading new model to GPU: " << path << std::endl;
        
        uint32_t currentTexOffset = static_cast<uint32_t>(m_globalTextures.size());
        MeshletModelBuilder builder;
        if (!builder.LoadAndSliceModel(path, currentTexOffset)) {
            throw std::runtime_error("Failed to load GLB: " + path);
        }
        
        // 1. 加载并追加纹理到全局无绑定 (Bindless) 数组中
        for (const auto& texData : builder.GetAllTextures()) {
            auto texture = std::make_unique<VulkanTexture>(m_device, m_resourceCommandPool.get(), texData.data(), texData.size());
            
            VkDescriptorImageInfo info{};
            info.imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
            info.imageView = texture->GetImageView();
            info.sampler = texture->GetSampler(); 
            
            m_globalImageInfos.push_back(info);
            m_globalTextures.push_back(std::move(texture));
        }

        VulkanModelResource res;
        res.totalMeshlets = static_cast<uint32_t>(builder.GetMeshlets().size());
        res.textureOffset = currentTexOffset;
        res.textureCount = static_cast<uint32_t>(builder.GetAllTextures().size());

        // 2. ★ 使用 RAII 方式创建 BDA Buffer，彻底告别裸指针管理！
        res.vertexBuffer   = CreateBDABuffer(builder.GetVertices());
        res.meshletBuffer  = CreateBDABuffer(builder.GetMeshlets());
        res.indexBuffer    = CreateBDABuffer(builder.GetMicroIndices());
        res.boundsBuffer   = CreateBDABuffer(builder.GetBounds());
        res.materialBuffer = CreateBDABuffer(builder.GetMaterials());

        // 缓存 BDA 物理地址，供 PushConstants 极速调用
        res.vAddr   = res.vertexBuffer ? res.vertexBuffer->GetDeviceAddress() : 0;
        res.mAddr   = res.meshletBuffer ? res.meshletBuffer->GetDeviceAddress() : 0;
        res.iAddr   = res.indexBuffer ? res.indexBuffer->GetDeviceAddress() : 0;
        res.bAddr   = res.boundsBuffer ? res.boundsBuffer->GetDeviceAddress() : 0;
        res.matAddr = res.materialBuffer ? res.materialBuffer->GetDeviceAddress() : 0;

        // 3. 构建供光追 (BLAS) 专用的全局索引
        const auto& vertices = builder.GetVertices();
        const auto& microIndices = builder.GetMicroIndices();
        const auto& meshlets = builder.GetMeshlets();

        std::vector<uint32_t> globalIndices;
        globalIndices.reserve(microIndices.size());
        for (const auto& m : meshlets) {
            for (uint32_t i = 0; i < m.triangleCount * 3; i++) {
                globalIndices.push_back(m.vertexOffset + microIndices[m.triangleOffset + i]); 
            }
        }

        uint32_t vertexCount = static_cast<uint32_t>(vertices.size());
        uint32_t indexCount = static_cast<uint32_t>(globalIndices.size());
        uint32_t vertexStride = vertices.empty() ? 0 : static_cast<uint32_t>(sizeof(vertices[0]));

        if (vertexCount > 0 && indexCount > 0) {
            res.globalIndexBuffer = CreateBDABuffer(globalIndices);
            res.globalIAddr = res.globalIndexBuffer->GetDeviceAddress();

            std::cout << "[VulkanBLAS] Triggering BLAS build for: " << path << std::endl;
            // 注意这里传入的是光追专用的全局索引 BDA (globalIAddr)
            res.blas = std::make_shared<VulkanBLAS>(
                m_device, 
                m_resourceCommandPool.get(), 
                res.vAddr, vertexCount, vertexStride,
                res.globalIAddr, indexCount 
            );
        }

        if (m_descriptorSet != VK_NULL_HANDLE && res.textureCount > 0) {
            VkWriteDescriptorSet write{};
            write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write.dstSet = m_descriptorSet;
            write.dstBinding = 0;
            // 从当前模型的新贴图偏移位置开始更新
            write.dstArrayElement = currentTexOffset; 
            write.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            write.descriptorCount = res.textureCount;
            write.pImageInfo = m_globalImageInfos.data() + currentTexOffset;
            vkUpdateDescriptorSets(m_device->GetNativeDevice(), 1, &write, 0, nullptr);
        }

        // 转移所有权到缓存中
        m_modelCache[path] = std::move(res);
        return m_modelCache[path];
    }


    static std::vector<char> ReadShaderFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::ate | std::ios::binary);
        if (!file.is_open()) throw std::runtime_error("Failed to open file: " + filename);
        size_t fileSize = (size_t)file.tellg();
        std::vector<char> buffer(fileSize);
        file.seekg(0); file.read(buffer.data(), fileSize); file.close();
        return buffer;
    }


    static VkShaderModule CreateShaderModule(VkDevice device, const std::vector<char>& code) {
        VkShaderModuleCreateInfo createInfo{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
        createInfo.codeSize = code.size();
        createInfo.pCode = reinterpret_cast<const uint32_t*>(code.data());
        VkShaderModule shaderModule;
        if (vkCreateShaderModule(device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) 
            throw std::runtime_error("Failed to create shader module!");
        return shaderModule;
    }


    void VulkanRenderingSystem::BuildPipelines() {
        std::cout << "[RenderingSystem] Building Pipelines..." << std::endl;
        const std::string SHADER_DIR = "C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/shader/";
        VkDevice device = m_device->GetNativeDevice();

        // =========================================================================
        // 阶段 1: 定义 Push Constants 范围
        // =========================================================================
        struct PushConstants {
            uint64_t frameDataAddr;    // 8 字节：包含所有全局相机、光照参数
            uint64_t modelDataAddr;    // 8 字节：包含当前物体的 Model 矩阵 (未来优化用)
            uint64_t vertexBuffer;     
            uint64_t meshletBuffer;    
            uint64_t indexBuffer;      
            uint64_t boundsBuffer;     
            uint64_t materialBuffer;   
            uint32_t totalMeshlets;    
            uint32_t textureOffset;
        };

        struct LightingPushConstants {
            uint64_t frameDataAddr; 
            uint64_t instanceDescAddr;   
            uint32_t stepSize;         // 4 字节：专门给 A-Trous 降噪用的步长
            uint32_t padding;          // 4 字节：补齐 8 字节对齐 (Vulkan 规范强烈建议)
        };

        VkPushConstantRange graphicsPushRange{}; 
        graphicsPushRange.stageFlags = VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT; 
        graphicsPushRange.size = sizeof(PushConstants); 

        VkPushConstantRange lightPushRange{}; 
        lightPushRange.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT; 
        lightPushRange.size = sizeof(LightingPushConstants);

        // =========================================================================
        // 阶段 2: 构建 Descriptor Sets (绑定数据源)
        // =========================================================================
        
        // 2.1 Bindless 贴图数组 (全局)
        VulkanDescriptorBuilder bindlessBuilder;
        VkDescriptorImageInfo dummyInfo = { m_gBufferSampler, m_gAlbedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };

        bindlessBuilder.BindImageArray(0, m_globalImageInfos.empty() ? &dummyInfo : m_globalImageInfos.data(), 1024, 
                                    m_globalImageInfos.empty() ? 1 : static_cast<uint32_t>(m_globalImageInfos.size()), 
                                    VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT, true);

        bindlessBuilder.Build(m_device, m_descriptorSetLayout, m_descriptorPool, m_descriptorSet);

        for (uint32_t ping = 0; ping < 2; ping++) {
            VkDevice dev = m_device->GetNativeDevice();

            // --- A. Lighting Pass (4个图片 + 1个TLAS) ---
            VkDescriptorSetLayoutBinding lightBindings[5] = {};
            for (uint32_t b = 0; b < 4; b++) {
                lightBindings[b].binding = b; lightBindings[b].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER; lightBindings[b].descriptorCount = 1; lightBindings[b].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }
            lightBindings[4].binding = 4; lightBindings[4].descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR; lightBindings[4].descriptorCount = 1; lightBindings[4].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            VkDescriptorSetLayoutCreateInfo lightLayoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            lightLayoutInfo.bindingCount = 5; lightLayoutInfo.pBindings = lightBindings;
            vkCreateDescriptorSetLayout(dev, &lightLayoutInfo, nullptr, &m_lightingLayouts[ping]);

            VkDescriptorPoolSize lightPoolSizes[2] = { {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 4}, {VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR, 1} };
            VkDescriptorPoolCreateInfo lightPoolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            lightPoolInfo.poolSizeCount = 2; lightPoolInfo.pPoolSizes = lightPoolSizes; lightPoolInfo.maxSets = 1;
            vkCreateDescriptorPool(dev, &lightPoolInfo, nullptr, &m_lightPools[ping]);

            VkDescriptorSetAllocateInfo lightAllocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            lightAllocInfo.descriptorPool = m_lightPools[ping]; lightAllocInfo.descriptorSetCount = 1; lightAllocInfo.pSetLayouts = &m_lightingLayouts[ping];
            vkAllocateDescriptorSets(dev, &lightAllocInfo, &m_lightingSets[ping]);

            // --- B. SVGF Temporal Pass (6个图片) ---
            VkDescriptorSetLayoutBinding temporalBindings[6] = {};
            for (uint32_t b = 0; b < 6; b++) {
                temporalBindings[b].binding = b; temporalBindings[b].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER; temporalBindings[b].descriptorCount = 1; temporalBindings[b].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }
            VkDescriptorSetLayoutCreateInfo temporalLayoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            temporalLayoutInfo.bindingCount = 6; temporalLayoutInfo.pBindings = temporalBindings;
            vkCreateDescriptorSetLayout(dev, &temporalLayoutInfo, nullptr, &m_svgfTemporalLayouts[ping]);
            
            VkDescriptorPoolSize tempPoolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 6};
            VkDescriptorPoolCreateInfo tempPoolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            tempPoolInfo.poolSizeCount = 1; tempPoolInfo.pPoolSizes = &tempPoolSize; tempPoolInfo.maxSets = 1;
            vkCreateDescriptorPool(dev, &tempPoolInfo, nullptr, &m_svgfTemporalPools[ping]);
            
            VkDescriptorSetAllocateInfo tempAllocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            tempAllocInfo.descriptorPool = m_svgfTemporalPools[ping]; tempAllocInfo.descriptorSetCount = 1; tempAllocInfo.pSetLayouts = &m_svgfTemporalLayouts[ping];
            vkAllocateDescriptorSets(dev, &tempAllocInfo, &m_svgfTemporalSets[ping]);

            // --- C. SVGF A-Trous Pass (4个图片，需要分4个Set) ---
            VkDescriptorSetLayoutBinding atrousBindings[4] = {};
            for (uint32_t b = 0; b < 4; b++) {
                atrousBindings[b].binding = b; atrousBindings[b].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER; atrousBindings[b].descriptorCount = 1; atrousBindings[b].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }
            VkDescriptorSetLayoutCreateInfo atrousLayoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            atrousLayoutInfo.bindingCount = 4; atrousLayoutInfo.pBindings = atrousBindings;
            vkCreateDescriptorSetLayout(dev, &atrousLayoutInfo, nullptr, &m_svgfATrousLayouts[ping]);
            
            VkDescriptorPoolSize atrousPoolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 16};
            VkDescriptorPoolCreateInfo atrousPoolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            atrousPoolInfo.poolSizeCount = 1; atrousPoolInfo.pPoolSizes = &atrousPoolSize; atrousPoolInfo.maxSets = 4;
            vkCreateDescriptorPool(dev, &atrousPoolInfo, nullptr, &m_svgfATrousPools[ping]);

            std::vector<VkDescriptorSetLayout> atrousLayoutsAlloc(4, m_svgfATrousLayouts[ping]);
            VkDescriptorSetAllocateInfo atrousAllocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            atrousAllocInfo.descriptorPool = m_svgfATrousPools[ping]; atrousAllocInfo.descriptorSetCount = 4; atrousAllocInfo.pSetLayouts = atrousLayoutsAlloc.data();
            vkAllocateDescriptorSets(dev, &atrousAllocInfo, m_svgfATrousSets[ping]);

            // --- D. TAA Pass (5个图片) ---
            VkDescriptorSetLayoutBinding taaBindings[5] = {};
            for (uint32_t b = 0; b < 5; b++) {
                taaBindings[b].binding = b; taaBindings[b].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER; taaBindings[b].descriptorCount = 1; taaBindings[b].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }
            VkDescriptorSetLayoutCreateInfo taaLayoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            taaLayoutInfo.bindingCount = 5; taaLayoutInfo.pBindings = taaBindings;
            vkCreateDescriptorSetLayout(dev, &taaLayoutInfo, nullptr, &m_taaLayouts[ping]);

            VkDescriptorPoolSize taaPoolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 5};
            VkDescriptorPoolCreateInfo taaPoolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            taaPoolInfo.poolSizeCount = 1; taaPoolInfo.pPoolSizes = &taaPoolSize; taaPoolInfo.maxSets = 1;
            vkCreateDescriptorPool(dev, &taaPoolInfo, nullptr, &m_taaPools[ping]);

            VkDescriptorSetAllocateInfo taaAllocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            taaAllocInfo.descriptorPool = m_taaPools[ping]; taaAllocInfo.descriptorSetCount = 1; taaAllocInfo.pSetLayouts = &m_taaLayouts[ping];
            vkAllocateDescriptorSets(dev, &taaAllocInfo, &m_taaSets[ping]);

            // --- E. Blit Pass (1个图片) ---
            VkDescriptorSetLayoutBinding blitBindings[1] = {};
            blitBindings[0].binding = 0; blitBindings[0].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER; blitBindings[0].descriptorCount = 1; blitBindings[0].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            VkDescriptorSetLayoutCreateInfo blitLayoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            blitLayoutInfo.bindingCount = 1; blitLayoutInfo.pBindings = blitBindings;
            vkCreateDescriptorSetLayout(dev, &blitLayoutInfo, nullptr, &m_blitLayouts[ping]);

            VkDescriptorPoolSize blitPoolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1};
            VkDescriptorPoolCreateInfo blitPoolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            blitPoolInfo.poolSizeCount = 1; blitPoolInfo.pPoolSizes = &blitPoolSize; blitPoolInfo.maxSets = 1;
            vkCreateDescriptorPool(dev, &blitPoolInfo, nullptr, &m_blitPools[ping]);

            VkDescriptorSetAllocateInfo blitAllocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            blitAllocInfo.descriptorPool = m_blitPools[ping]; blitAllocInfo.descriptorSetCount = 1; blitAllocInfo.pSetLayouts = &m_blitLayouts[ping];
            vkAllocateDescriptorSets(dev, &blitAllocInfo, &m_blitSets[ping]);
        }

        UpdateDescriptorSets();

        // =========================================================================
        // 阶段 3: 加载 Shader Modules
        // =========================================================================
        // 3.1 Mesh Shaders
        VkShaderModule taskShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "car_task.spv"));
        VkShaderModule meshShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "car_mesh.spv"));
        VkShaderModule fragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "car_frag.spv"));

        // 3.2 Fullscreen Quad 相关的 Shaders (Lighting, Denoise, TAA, Blit 共用 Vert)
        VkShaderModule fsVertShader    = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "lighting_vert.spv"));
        VkShaderModule lightFragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "lighting_frag.spv"));
        VkShaderModule svgfTemporalFragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "svgf_temporal_frag.spv"));
        VkShaderModule svgfATrousFragShader   = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "svgf_atrous_frag.spv"));
        VkShaderModule taaFragShader   = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "taa_frag.spv"));
        VkShaderModule blitFragShader  = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "blit_frag.spv"));

        // =========================================================================
        // 阶段 4: 构建 Pipeline Layouts & Pipelines
        // =========================================================================
        
        // 4.1 Graphics Pipeline (G-Buffer)
        VkPipelineLayoutCreateInfo gPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        gPipeLayoutInfo.pushConstantRangeCount = 1; 
        gPipeLayoutInfo.pPushConstantRanges = &graphicsPushRange; 
        gPipeLayoutInfo.setLayoutCount = 1; 
        gPipeLayoutInfo.pSetLayouts = &m_descriptorSetLayout;
        vkCreatePipelineLayout(device, &gPipeLayoutInfo, nullptr, &m_graphicsPipelineLayout);

        m_graphicsPipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_TASK_BIT_EXT, taskShader)
            .AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(true, true, VK_COMPARE_OP_GREATER_OR_EQUAL)
            .DisableColorBlendAttachments(3) // Albedo, Normal, Velocity
            .SetPipelineLayout(m_graphicsPipelineLayout)
            .Build(m_device, { VK_FORMAT_R8G8B8A8_UNORM, VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16_SFLOAT }, VK_FORMAT_D32_SFLOAT); 

        // 4.2 Lighting Pipeline
        VkDescriptorSetLayout lightSetLayouts[2] = { m_lightingLayouts[0], m_descriptorSetLayout }; // 需要 Bindless 材质支持
        VkPipelineLayoutCreateInfo lightPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        lightPipeLayoutInfo.setLayoutCount = 2; 
        lightPipeLayoutInfo.pSetLayouts = lightSetLayouts; 
        lightPipeLayoutInfo.pushConstantRangeCount = 1; 
        lightPipeLayoutInfo.pPushConstantRanges = &lightPushRange; 
        vkCreatePipelineLayout(device, &lightPipeLayoutInfo, nullptr, &m_lightingPipelineLayout);

        m_lightingPipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, fsVertShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, lightFragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false).AddColorBlendAttachment(false) // DirectLight, NoisyGI
            .SetPipelineLayout(m_lightingPipelineLayout)
            .Build(m_device, { VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED);

        // 4.3 SVGF Temporal Pipeline (输出 2 个附件：GI History 和 Moments)
        VkPipelineLayoutCreateInfo temporalPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        temporalPipeLayoutInfo.setLayoutCount = 1; 
        temporalPipeLayoutInfo.pSetLayouts = &m_svgfTemporalLayouts[0]; 
        temporalPipeLayoutInfo.pushConstantRangeCount = 1; 
        temporalPipeLayoutInfo.pPushConstantRanges = &lightPushRange; 
        vkCreatePipelineLayout(device, &temporalPipeLayoutInfo, nullptr, &m_svgfTemporalPipelineLayout);

        m_svgfTemporalPipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, fsVertShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, svgfTemporalFragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            // ★ 注意：这里必须 Add 两次，因为输出给 m_gGIHistory 和 m_gMomentsHistory 两个附件
            .AddColorBlendAttachment(false) 
            .AddColorBlendAttachment(false) 
            .SetPipelineLayout(m_svgfTemporalPipelineLayout)
            .Build(m_device, { VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16_SFLOAT }, VK_FORMAT_UNDEFINED);

        // 4.4 SVGF A-Trous Pipeline (输出 1 个附件：DenoisedGI)
        VkPipelineLayoutCreateInfo atrousPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        atrousPipeLayoutInfo.setLayoutCount = 1; 
        atrousPipeLayoutInfo.pSetLayouts = &m_svgfATrousLayouts[0]; 
        atrousPipeLayoutInfo.pushConstantRangeCount = 1; 
        atrousPipeLayoutInfo.pPushConstantRanges = &lightPushRange; 
        vkCreatePipelineLayout(device, &atrousPipeLayoutInfo, nullptr, &m_svgfATrousPipelineLayout);

        m_svgfATrousPipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, fsVertShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, svgfATrousFragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false) // 只输出 DenoisedGI
            .SetPipelineLayout(m_svgfATrousPipelineLayout)
            .Build(m_device, { VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED);

        // 4.5 TAA Pipeline
        VkPipelineLayoutCreateInfo taaPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        taaPipeLayoutInfo.setLayoutCount = 1; 
        taaPipeLayoutInfo.pSetLayouts = &m_taaLayouts[0]; 
        taaPipeLayoutInfo.pushConstantRangeCount = 1;               // [+] 声明 Push Constant 数量
        taaPipeLayoutInfo.pPushConstantRanges = &lightPushRange;    // [+] 绑定 Push Constant 范围
        vkCreatePipelineLayout(device, &taaPipeLayoutInfo, nullptr, &m_taaPipelineLayout);

        m_taaPipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, fsVertShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, taaFragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false)
            .SetPipelineLayout(m_taaPipelineLayout)
            .Build(m_device, { VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED);

        // 4.6 Blit Pipeline (推向屏幕)
        VkPipelineLayoutCreateInfo blitPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        blitPipeLayoutInfo.setLayoutCount = 1; 
        blitPipeLayoutInfo.pSetLayouts = &m_blitLayouts[0]; 
        blitPipeLayoutInfo.pushConstantRangeCount = 1;              // [+] 声明 Push Constant 数量
        blitPipeLayoutInfo.pPushConstantRanges = &lightPushRange;   // [+] 绑定 Push Constant 范围
        vkCreatePipelineLayout(device, &blitPipeLayoutInfo, nullptr, &m_blitPipelineLayout);

        m_blitPipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, fsVertShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, blitFragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false)
            .SetPipelineLayout(m_blitPipelineLayout)
            .Build(m_device, { m_renderer->GetSwapchainFormat() }, VK_FORMAT_D32_SFLOAT); // 目标是屏幕

        // =========================================================================
        // 阶段 5: 清理 Shader Modules
        // =========================================================================
        vkDestroyShaderModule(device, taskShader, nullptr);
        vkDestroyShaderModule(device, meshShader, nullptr);
        vkDestroyShaderModule(device, fragShader, nullptr);
        vkDestroyShaderModule(device, fsVertShader, nullptr);
        vkDestroyShaderModule(device, lightFragShader, nullptr);
        vkDestroyShaderModule(device, svgfTemporalFragShader, nullptr);
        vkDestroyShaderModule(device, svgfATrousFragShader, nullptr);
        vkDestroyShaderModule(device, taaFragShader, nullptr);
        vkDestroyShaderModule(device, blitFragShader, nullptr);

        std::cout << "[RenderingSystem] Pipelines Built Successfully." << std::endl;
    }


    void VulkanRenderingSystem::Tick(Registry& registry, float deltaTime) {
        m_frameIndex++;

        if (!m_CmdDrawMeshTasksEXT) {
            m_CmdDrawMeshTasksEXT = (PFN_vkCmdDrawMeshTasksEXT)vkGetDeviceProcAddr(m_device->GetNativeDevice(), "vkCmdDrawMeshTasksEXT");
        }

        // ====================================================================
        // ★ 回归纯粹：最完美的 TAA Jitter 计算 (直接基于物理分辨率 m_width/m_height)
        // ====================================================================
        uint32_t jitterIndex = m_frameIndex % 8 + 1; 
        float jitterX = CreateHaltonSequence(jitterIndex, 2) - 0.5f; 
        float jitterY = CreateHaltonSequence(jitterIndex, 3) - 0.5f;
        float ndcJitterX = (jitterX * 2.0f) / static_cast<float>(m_width > 0 ? m_width : 1); 
        float ndcJitterY = (jitterY * 2.0f) / static_cast<float>(m_height > 0 ? m_height : 1);

        uint32_t ping = m_frameIndex % 2;       
        VkCommandBuffer cmd = m_renderer->BeginFrame();
        if (cmd == VK_NULL_HANDLE) return; 

        // ★ 回归纯粹：唯一的全屏视口和裁剪框
        VkViewport viewport{ 0.0f, 0.0f, (float)m_width, (float)m_height, 0.0f, 1.0f };
        VkRect2D scissor{ {0, 0}, {m_width, m_height} };

        Matrix4x4 viewMat, projMatUnjittered, projMatJittered;
        Vector3 cameraPos;
        auto cameraView = registry.view<TransformComponent, CameraComponent>();
        
        for (auto entity : cameraView) {
            auto& camera = cameraView.get<CameraComponent>(entity);
            auto& transform = cameraView.get<TransformComponent>(entity);
            cameraPos = transform.getPosition();

            // 直接使用真实窗口宽高比
            float aspect = static_cast<float>(m_width) / static_cast<float>(m_height > 0 ? m_height : 1);
            Matrix4x4 baseProj = camera.BuildPerspectiveInfiniteReverseZ(camera.getFov(), aspect, camera.getzNear());
            camera.setProjectionMatrix(baseProj);

            viewMat = camera.getViewMatrix(); 
            projMatUnjittered = camera.getProjectionMatrix();
            projMatUnjittered[1][1] *= -1.0f; // Vulkan Y-flip
            
            projMatJittered = projMatUnjittered;
            projMatJittered[0][2] += ndcJitterX; 
            projMatJittered[1][2] += ndcJitterY;
            break; 
        }

        Lizeral::Vector3 lightDir(1.0f, 0.5f, 1.0f);
        Lizeral::Vector3 lightColor(1.0f, 0.85f, 0.7f);
        float lightIntensity = 4.0f;

        auto lightView = registry.view<TransformComponent, DirectionLightComponent>();
        for (auto entity : lightView) {
            auto& trans = lightView.get<TransformComponent>(entity);
            auto& light = lightView.get<DirectionLightComponent>(entity);
            lightDir = trans.getForward().normalize();
            lightColor = light.getColor(); 
            lightIntensity = light.getIntensity();
            break; // 通常只有一个主光源
        }

        Matrix4x4 currentVP = projMatUnjittered * viewMat;

        GlobalFrameData frameData{};
        // ★ 修复：将所有相机矩阵转置 (transpose) 后再喂给 VRAM！
        frameData.viewProj = currentVP.transpose();
        frameData.invViewProj = currentVP.inverse().transpose();
        frameData.prevViewProj = m_prevViewProj.transpose();
        
        frameData.cameraPos = cameraPos;
        frameData.lightDir = lightDir;
        frameData.lightColor = lightColor;
        frameData.lightIntensity = lightIntensity; 
        frameData.frameIndex = m_frameIndex;
        frameData.jitter = Vector2(ndcJitterX, ndcJitterY);

        if (!m_frameResource.buffer) {
            std::vector<GlobalFrameData> dummy(1);
            m_frameResource.buffer = CreateBDABuffer(dummy);
            m_frameResource.addr = m_frameResource.buffer->GetDeviceAddress();
        }

        m_frameResource.buffer->WriteData(&frameData, sizeof(GlobalFrameData));

        // ====================================================================
        // 步骤 A: 收集实例，构建动态 TLAS 光追树
        // ====================================================================
        struct RTInstanceDesc {
            uint64_t vertexBuffer; uint64_t indexBuffer; uint64_t materialBuffer;
            uint32_t textureOffset; uint32_t padding[3];
        };
        RTInstanceDesc* mappedDesc = static_cast<RTInstanceDesc*>(m_rtInstanceBuffer->Map());

        std::vector<VkAccelerationStructureInstanceKHR> tlasInstances;
        uint32_t customInstanceId = 0;

        auto modelView = registry.view<TransformComponent, VulkanModelComponent>();
        for (auto entity : modelView) {
            auto& transform = modelView.get<TransformComponent>(entity);
            auto& modelComp = modelView.get<VulkanModelComponent>(entity);
            
            if (modelComp.getVulkanModelPath().empty()) continue;

            auto& res = GetOrLoadModel(modelComp.getVulkanModelPath());
            if (!res.IsValid() || !res.blas) continue;

            mappedDesc[customInstanceId].vertexBuffer = res.vAddr;
            mappedDesc[customInstanceId].indexBuffer  = res.globalIAddr;
            mappedDesc[customInstanceId].materialBuffer = res.matAddr;
            mappedDesc[customInstanceId].textureOffset  = res.textureOffset;

            Matrix4x4 modelMat = transform.getMatrix(); 
            VkTransformMatrixKHR vkTransform{};
            vkTransform.matrix[0][0] = modelMat[0][0]; vkTransform.matrix[0][1] = modelMat[0][1]; vkTransform.matrix[0][2] = modelMat[0][2]; vkTransform.matrix[0][3] = modelMat[0][3];
            vkTransform.matrix[1][0] = modelMat[1][0]; vkTransform.matrix[1][1] = modelMat[1][1]; vkTransform.matrix[1][2] = modelMat[1][2]; vkTransform.matrix[1][3] = modelMat[1][3];
            vkTransform.matrix[2][0] = modelMat[2][0]; vkTransform.matrix[2][1] = modelMat[2][1]; vkTransform.matrix[2][2] = modelMat[2][2]; vkTransform.matrix[2][3] = modelMat[2][3];

            VkAccelerationStructureInstanceKHR instance{};
            instance.transform = vkTransform; 
            instance.instanceCustomIndex = customInstanceId++; 
            instance.mask = 0xFF; 
            instance.flags = VK_GEOMETRY_INSTANCE_TRIANGLE_FACING_CULL_DISABLE_BIT_KHR; 
            instance.accelerationStructureReference = res.blas->GetDeviceAddress();
            tlasInstances.push_back(instance);
        }
        m_rtInstanceBuffer->Unmap();

        if (!tlasInstances.empty()) {
            m_tlas->Build(cmd, ping, tlasInstances, false); // 根据你的 TLAS 更新逻辑可能需要调整 false
            VkMemoryBarrier memoryBarrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
            memoryBarrier.srcAccessMask = VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR;
            memoryBarrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
            vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
        }

        // ====================================================================
        // 步骤 B: 极速热更新 Descriptor
        // ====================================================================
        VkAccelerationStructureKHR currentTlas = m_tlas->GetHandle(ping);
        if (currentTlas != VK_NULL_HANDLE) {
            VkWriteDescriptorSetAccelerationStructureKHR asWrite{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET_ACCELERATION_STRUCTURE_KHR};
            asWrite.accelerationStructureCount = 1; 
            asWrite.pAccelerationStructures = &currentTlas;

            VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
            write.dstSet = m_lightingSets[ping]; 
            write.dstBinding = 4; 
            write.dstArrayElement = 0;
            write.descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
            write.descriptorCount = 1; 
            write.pNext = &asWrite;
            
            vkUpdateDescriptorSets(m_device->GetNativeDevice(), 1, &write, 0, nullptr);
        }

        // ====================================================================
        // 步骤 C: G-Buffer Pass
        // ====================================================================
        VkImageLayout currentLayout = m_firstFrame ? VK_IMAGE_LAYOUT_UNDEFINED : VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        TransitionImageLayout(cmd, m_gAlbedoMetallic.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        TransitionImageLayout(cmd, m_gNormalRoughness.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        TransitionImageLayout(cmd, m_gDepth.image, currentLayout, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);
        TransitionImageLayout(cmd, m_gVelocity.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

        VkRenderingAttachmentInfo colorAttachments[3] = {};
        colorAttachments[0].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO; colorAttachments[0].imageView = m_gAlbedoMetallic.view; colorAttachments[0].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; colorAttachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; colorAttachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE; colorAttachments[0].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};
        colorAttachments[1].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO; colorAttachments[1].imageView = m_gNormalRoughness.view; colorAttachments[1].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; colorAttachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; colorAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_STORE; colorAttachments[1].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};
        colorAttachments[2].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO; colorAttachments[2].imageView = m_gVelocity.view; colorAttachments[2].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; colorAttachments[2].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; colorAttachments[2].storeOp = VK_ATTACHMENT_STORE_OP_STORE; colorAttachments[2].clearValue.color = {0.0f, 0.0f, 0.0f, 0.0f};
        VkRenderingAttachmentInfo depthAttachment{VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO}; depthAttachment.imageView = m_gDepth.view; depthAttachment.imageLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL; depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR; depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE; depthAttachment.clearValue.depthStencil = {0.0f, 0};

        VkRenderingInfo renderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO}; 
        renderInfo.renderArea = scissor; renderInfo.layerCount = 1; 
        renderInfo.colorAttachmentCount = 3; renderInfo.pColorAttachments = colorAttachments; renderInfo.pDepthAttachment = &depthAttachment;
        
        vkCmdBeginRendering(cmd, &renderInfo);
        vkCmdSetViewport(cmd, 0, 1, &viewport); vkCmdSetScissor(cmd, 0, 1, &scissor);
        
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_graphicsPipeline);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_graphicsPipelineLayout, 0, 1, &m_descriptorSet, 0, nullptr);

        Matrix4x4 currentViewProjJittered = projMatJittered * viewMat;
        Matrix4x4 currentViewProjUnjittered = projMatUnjittered * viewMat;
        if (m_isFirstFrameRun) m_prevViewProj = currentViewProjUnjittered;

        struct PushConstants {
            uint64_t frameDataAddr;    // 8 字节：包含所有全局相机、光照参数
            uint64_t modelDataAddr;    // 8 字节：包含当前物体的 Model 矩阵 (未来优化用)
            uint64_t vertexBuffer;     
            uint64_t meshletBuffer;    
            uint64_t indexBuffer;      
            uint64_t boundsBuffer;     
            uint64_t materialBuffer;   
            uint32_t totalMeshlets;    
            uint32_t textureOffset;
        };

        GPUInstanceData* mappedInstanceData = static_cast<GPUInstanceData*>(m_globalInstanceBuffer->Map());
        uint64_t baseInstanceAddr = m_globalInstanceBuffer->GetDeviceAddress();
        uint32_t currentInstanceIndex = 0;

        for (auto entity : modelView) {
            auto& transform = modelView.get<TransformComponent>(entity);
            auto& modelComp = modelView.get<VulkanModelComponent>(entity);
            if (modelComp.getVulkanModelPath().empty()) continue;

            const auto& res = GetOrLoadModel(modelComp.getVulkanModelPath());
            if (!res.IsValid()) continue;

            Matrix4x4 currentModel = transform.getMatrix();
            uint32_t entityId = static_cast<uint32_t>(entity);
            Matrix4x4 prevModel = m_isFirstFrameRun ? currentModel : m_prevModelMats[entityId];

            PushConstants pushData{};

            if (currentInstanceIndex >= 10000) {
                std::cerr << "WARNING: Max instances (10000) reached!" << std::endl;
                break;
            }

            // 2. 疯狂写入数据
            mappedInstanceData[currentInstanceIndex].Model = currentModel.transpose();
            mappedInstanceData[currentInstanceIndex].prevModel = prevModel.transpose();
            
            pushData.frameDataAddr = m_frameResource.addr;
            pushData.modelDataAddr = baseInstanceAddr + (currentInstanceIndex * sizeof(GPUInstanceData));
            pushData.vertexBuffer = res.vAddr; pushData.meshletBuffer = res.mAddr; pushData.indexBuffer = res.iAddr; 
            pushData.boundsBuffer = res.bAddr; pushData.materialBuffer = res.matAddr;
            pushData.totalMeshlets = res.totalMeshlets; pushData.textureOffset = res.textureOffset; 

            vkCmdPushConstants(cmd, m_graphicsPipelineLayout, VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(PushConstants), &pushData);
            m_prevModelMats[entityId] = currentModel;
            m_CmdDrawMeshTasksEXT(cmd, (res.totalMeshlets + 63) / 64, 1, 1);
            currentInstanceIndex++;
        }
        vkCmdEndRendering(cmd);

        m_prevViewProj = currentViewProjUnjittered;

        // ====================================================================
        // 步骤 D: 全屏后处理链
        // ====================================================================
        if (m_isFirstFrameRun) {
            TransitionImageLayout(cmd, m_gHistory[0].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_gHistory[1].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_gGIHistory[0].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_gGIHistory[1].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_gMomentsHistory[0].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_gMomentsHistory[1].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        }
        
        TransitionImageLayout(cmd, m_gAlbedoMetallic.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        TransitionImageLayout(cmd, m_gNormalRoughness.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        TransitionImageLayout(cmd, m_gDepth.image, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);
        TransitionImageLayout(cmd, m_gVelocity.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

        struct LightingPushConstants {
            uint64_t frameDataAddr; 
            uint64_t instanceDescAddr;   
            uint32_t stepSize;         // 4 字节：专门给 A-Trous 降噪用的步长
            uint32_t padding;          // 4 字节：补齐 8 字节对齐 (Vulkan 规范强烈建议)
        };

        LightingPushConstants lightPc{};
        lightPc.frameDataAddr = m_frameResource.addr;
        lightPc.instanceDescAddr = m_rtInstanceBuffer->GetDeviceAddress();


        auto DrawFullscreenPass = [&](VkPipeline pipeline, VkPipelineLayout layout, VkDescriptorSet set, GBufferAttachment* outputs, uint32_t outputCount, bool bindExtraSet = false) {
            VkRenderingAttachmentInfo attInfos[2] = {};
            for(uint32_t i=0; i<outputCount; i++) {
                TransitionImageLayout(cmd, outputs[i].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
                attInfos[i].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
                attInfos[i].imageView = outputs[i].view;
                attInfos[i].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                attInfos[i].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
                attInfos[i].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
            }
            VkRenderingInfo rInfo{VK_STRUCTURE_TYPE_RENDERING_INFO}; 
            rInfo.renderArea = scissor; rInfo.layerCount = 1; rInfo.colorAttachmentCount = outputCount; rInfo.pColorAttachments = attInfos;
            
            vkCmdBeginRendering(cmd, &rInfo);
            vkCmdSetViewport(cmd, 0, 1, &viewport); vkCmdSetScissor(cmd, 0, 1, &scissor);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline);
            if (bindExtraSet) {
                VkDescriptorSet boundSets[2] = { set, m_descriptorSet }; 
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, layout, 0, 2, boundSets, 0, nullptr);
            } else {
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, layout, 0, 1, &set, 0, nullptr);
            }
            vkCmdPushConstants(cmd, layout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(LightingPushConstants), &lightPc);
            vkCmdDraw(cmd, 3, 1, 0, 0);
            vkCmdEndRendering(cmd);
            
            for(uint32_t i=0; i<outputCount; i++) {
                TransitionImageLayout(cmd, outputs[i].image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            }
        };

        // 1. Lighting Pass
        GBufferAttachment lightOuts[2] = { m_gDirectLight, m_gNoisyGI };
        DrawFullscreenPass(m_lightingPipeline, m_lightingPipelineLayout, m_lightingSets[ping], lightOuts, 2, true);

        // 2. SVGF Temporal Pass
        GBufferAttachment temporalOuts[2] = { m_gGIHistory[ping], m_gMomentsHistory[ping] };
        DrawFullscreenPass(m_svgfTemporalPipeline, m_svgfTemporalPipelineLayout, m_svgfTemporalSets[ping], temporalOuts, 2);

        // 3. SVGF A-Trous Pass
        if (m_isFirstFrameRun) {
            TransitionImageLayout(cmd, m_gDenoisedGITemp.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        }

        int stepSizes[] = { 1, 2, 4, 8 };
        GBufferAttachment* currentOutput = &m_gDenoisedGITemp;

        for (int i = 0; i < 4; i++) {
            lightPc.stepSize = stepSizes[i];
            DrawFullscreenPass(m_svgfATrousPipeline, m_svgfATrousPipelineLayout, m_svgfATrousSets[ping][i], currentOutput, 1);
            currentOutput = (currentOutput == &m_gDenoisedGITemp) ? &m_gDenoisedGI : &m_gDenoisedGITemp;
        }

        // 4. TAA Pass
        DrawFullscreenPass(m_taaPipeline, m_taaPipelineLayout, m_taaSets[ping], &m_gHistory[ping], 1);

        // 5. Blit Pass
        m_renderer->BeginRendering(cmd); 
        vkCmdSetViewport(cmd, 0, 1, &viewport); vkCmdSetScissor(cmd, 0, 1, &scissor);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_blitPipeline);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_blitPipelineLayout, 0, 1, &m_blitSets[ping], 0, nullptr);
        vkCmdDraw(cmd, 3, 1, 0, 0);
        m_renderer->EndRendering(cmd);

        m_renderer->EndFrame();
        m_isFirstFrameRun = m_firstFrame = false;
    }


    void VulkanRenderingSystem::UpdateDescriptorSets() {
        VkDevice device = m_device->GetNativeDevice();

        VkDescriptorImageInfo gInfos[4] = {
            { m_gBufferSampler, m_gAlbedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { m_gBufferSampler, m_gNormalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { m_gBufferSampler, m_gDepth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { m_gBufferSampler, m_gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
        };

        for (uint32_t i = 0; i < 2; i++) {
            uint32_t ping = i;
            uint32_t pong = (i + 1) % 2;

            // 1. 更新 Lighting Pass Sets (前4个槽位)
            VkWriteDescriptorSet lightWrites[4] = {};
            for (uint32_t b = 0; b < 4; b++) {
                lightWrites[b].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                lightWrites[b].dstSet = m_lightingSets[ping];
                lightWrites[b].dstBinding = b;
                lightWrites[b].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                lightWrites[b].descriptorCount = 1;
                lightWrites[b].pImageInfo = &gInfos[b];
            }
            vkUpdateDescriptorSets(device, 4, lightWrites, 0, nullptr);

            // 2. 更新 SVGF Temporal Pass Sets
            VkDescriptorImageInfo temporalInfos[6] = {
                { m_gBufferSampler, m_gNoisyGI.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { m_gBufferSampler, m_gNormalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { m_gBufferSampler, m_gDepth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { m_gBufferSampler, m_gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { m_gBufferSampler, m_gGIHistory[pong].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { m_gBufferSampler, m_gMomentsHistory[pong].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
            };
            VkWriteDescriptorSet temporalWrites[6] = {};
            for (uint32_t b = 0; b < 6; b++) {
                temporalWrites[b].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                temporalWrites[b].dstSet = m_svgfTemporalSets[ping];
                temporalWrites[b].dstBinding = b;
                temporalWrites[b].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                temporalWrites[b].descriptorCount = 1;
                temporalWrites[b].pImageInfo = &temporalInfos[b];
            }
            vkUpdateDescriptorSets(device, 6, temporalWrites, 0, nullptr);

            // 3. 更新 SVGF A-Trous Pass Sets (4轮迭代)
            VkImageView atrousInputs[4] = { m_gGIHistory[ping].view, m_gDenoisedGITemp.view, m_gDenoisedGI.view, m_gDenoisedGITemp.view };
            for (int iter = 0; iter < 4; iter++) {
                VkDescriptorImageInfo atrousInfos[4] = {
                    { m_gBufferSampler, atrousInputs[iter], VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                    { m_gBufferSampler, m_gMomentsHistory[ping].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                    { m_gBufferSampler, m_gNormalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                    { m_gBufferSampler, m_gDepth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
                };
                VkWriteDescriptorSet atrousWrites[4] = {};
                for (uint32_t b = 0; b < 4; b++) {
                    atrousWrites[b].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                    atrousWrites[b].dstSet = m_svgfATrousSets[ping][iter];
                    atrousWrites[b].dstBinding = b;
                    atrousWrites[b].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                    atrousWrites[b].descriptorCount = 1;
                    atrousWrites[b].pImageInfo = &atrousInfos[b];
                }
                vkUpdateDescriptorSets(device, 4, atrousWrites, 0, nullptr);
            }

            // 4. 更新 TAA Pass Sets
            VkDescriptorImageInfo taaInfos[5] = { 
                { m_gBufferSampler, m_gDenoisedGI.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { m_gBufferSampler, m_gHistory[pong].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }, 
                { m_gBufferSampler, m_gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { m_gBufferSampler, m_gAlbedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }, 
                { m_gBufferSampler, m_gDirectLight.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }     
            };
            VkWriteDescriptorSet taaWrites[5] = {};
            for (uint32_t b = 0; b < 5; b++) {
                taaWrites[b].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                taaWrites[b].dstSet = m_taaSets[ping];
                taaWrites[b].dstBinding = b;
                taaWrites[b].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                taaWrites[b].descriptorCount = 1;
                taaWrites[b].pImageInfo = &taaInfos[b];
            }
            vkUpdateDescriptorSets(device, 5, taaWrites, 0, nullptr);

            // 5. 更新 Blit Pass Sets
            VkDescriptorImageInfo blitInfos[1] = { { m_gBufferSampler, m_gHistory[ping].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL } };
            VkWriteDescriptorSet blitWrite{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
            blitWrite.dstSet = m_blitSets[ping];
            blitWrite.dstBinding = 0;
            blitWrite.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            blitWrite.descriptorCount = 1;
            blitWrite.pImageInfo = &blitInfos[0];
            vkUpdateDescriptorSets(device, 1, &blitWrite, 0, nullptr);
        }
    }

} // namespace Lizeral