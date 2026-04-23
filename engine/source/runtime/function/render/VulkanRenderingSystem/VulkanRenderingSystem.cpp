#include "VulkanRenderingSystem.h"
#include "runtime/function/render/rhi/vulkan/VulkanPipelineBuilder.h"
#include "runtime/function/render/rhi/vulkan/VulkanDescriptorBuilder.h"
#include "runtime/function/render/MeshletBuilder/MeshletModelBuilder.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Model/VulkanModelComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Light/DirectionalLightComponent.h"
#include "runtime/function/ecs/components/Light/PointLightComponent.h"

#include <iostream>
#include <fstream>
#include <algorithm>

namespace Lizeral {

    namespace {
        struct EditorGridPushConstants {
            Matrix4x4 invViewProj;
            Vector4 cameraPosAndPlaneHeight;
            Vector4 viewportSizeAndSpacing;
            Vector4 fadeAndOpacity;
        };
    } // namespace

    VulkanRenderingSystem::LightingProfile VulkanRenderingSystem::ResolveLightingProfile(RenderPipelinePreset preset) const {
        switch (preset) {
            case RenderPipelinePreset::SSGI:
                return LightingProfile{1, 0};
            case RenderPipelinePreset::RTGI:
                return LightingProfile{2, 1};
            case RenderPipelinePreset::Stable:
            default:
                return LightingProfile{0, 0};
        }
    }

    void VulkanRenderingSystem::SetRenderPipelinePreset(RenderPipelinePreset preset) {
        if (m_renderPreset == preset) return;

        m_renderPreset = preset;
        if (!m_useManualLightingProfile) {
            m_lightingProfile = ResolveLightingProfile(preset);
        }

        // Allow preset switch at runtime by rebuilding descriptor state and pipelines.
        if (m_device != nullptr) {
            VkDevice device = m_device->GetNativeDevice();
            vkDeviceWaitIdle(device);
            DestroyPipelines(device);
            DestroyDescriptors(device);
            BuildPipelines();
            InvalidateTemporalHistory();
        }
    }

    void VulkanRenderingSystem::SetLightingProfile(const LightingProfile& profile) {
        const bool unchanged =
            m_useManualLightingProfile &&
            m_lightingProfile.giQualityLevel == profile.giQualityLevel &&
            m_lightingProfile.shadowQualityLevel == profile.shadowQualityLevel;
        if (unchanged) return;

        m_useManualLightingProfile = true;
        m_lightingProfile = profile;

        if (m_device != nullptr) {
            VkDevice device = m_device->GetNativeDevice();
            vkDeviceWaitIdle(device);
            DestroyPipelines(device);
            DestroyDescriptors(device);
            BuildPipelines();
            InvalidateTemporalHistory();
        }
    }

    void VulkanRenderingSystem::ResetLightingProfileToPreset() {
        m_useManualLightingProfile = false;
        SetLightingProfile(ResolveLightingProfile(m_renderPreset));
        m_useManualLightingProfile = false;
    }

    template <typename T>
    std::unique_ptr<VulkanBuffer> VulkanRenderingSystem::CreateBDABuffer(const std::vector<T>& data) {
        if (data.empty()) return nullptr;

        VkDeviceSize bufferSize = data.size() * sizeof(T);
        
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
        if (!m_useManualLightingProfile) {
            m_lightingProfile = ResolveLightingProfile(m_renderPreset);
        }

        std::cout << "[VulkanRenderingSystem] Initializing..." << std::endl;
        std::cout
            << "[VulkanRenderingSystem] Lighting profile: GI=" << m_lightingProfile.giQualityLevel
            << " Shadow=" << m_lightingProfile.shadowQualityLevel
            << (m_useManualLightingProfile ? " (manual)" : " (preset)")
            << std::endl;

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

            std::vector<GPUPointLight> dummyLights(1024);
            m_pointLightBuffer = CreateBDABuffer(dummyLights);
        }

        CreateAttachments();

        VkSamplerCreateInfo samplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
        samplerInfo.magFilter = VK_FILTER_LINEAR; samplerInfo.minFilter = VK_FILTER_LINEAR;
        samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
        samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE; 
        samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE; 
        samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
        vkCreateSampler(m_device->GetNativeDevice(), &samplerInfo, nullptr, &m_gBufferSampler);

        m_debugLineBuffer = std::make_unique<VulkanBuffer>(
            m_device, m_maxDebugLineBufferSize, 
            VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, 
            VMA_MEMORY_USAGE_CPU_TO_GPU
        );

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
        vkDeviceWaitIdle(m_device->GetNativeDevice());
        VkDevice device = m_device->GetNativeDevice();

        std::cout << "[RenderingSystem] Shutting down..." << std::endl;

        m_modelCache.clear();      
        m_globalTextures.clear();  
        m_globalImageInfos.clear();
        m_globalTexturePathCache.clear();
        m_overrideMaterialBuffers.clear();

        if (m_rtInstanceBuffer) m_rtInstanceBuffer.reset();
        if (m_globalInstanceBuffer) m_globalInstanceBuffer.reset(); 
        if (m_frameResource.buffer) m_frameResource.buffer.reset(); 
        if (m_pointLightBuffer) m_pointLightBuffer.reset();
        if (m_debugLineBuffer) m_debugLineBuffer.reset();

        DestroyPipelines(device);
        DestroyDescriptors(device);
        DestroyAttachments();

        if (m_gBufferSampler) vkDestroySampler(device, m_gBufferSampler, nullptr);
        if (m_tlas) m_tlas.reset();
        if (m_resourceCommandPool) m_resourceCommandPool.reset();

        std::cout << "[RenderingSystem] Shutdown complete." << std::endl;
    }

    void VulkanRenderingSystem::DestroyPipelines(VkDevice device){
        if (m_blitPipeline)       vkDestroyPipeline(device, m_blitPipeline, nullptr);
        if (m_blitPipelineLayout) vkDestroyPipelineLayout(device, m_blitPipelineLayout, nullptr);

        if (m_transparentPipeline)       vkDestroyPipeline(device, m_transparentPipeline, nullptr);
        if (m_transparentPipelineLayout) vkDestroyPipelineLayout(device, m_transparentPipelineLayout, nullptr);

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

        if (m_editorGridPipeline) vkDestroyPipeline(device, m_editorGridPipeline, nullptr);
        if (m_editorGridPipelineLayout) vkDestroyPipelineLayout(device, m_editorGridPipelineLayout, nullptr);

        if (m_debugLinePipeline) vkDestroyPipeline(device, m_debugLinePipeline, nullptr);
        if (m_debugLinePipelineLayout) vkDestroyPipelineLayout(device, m_debugLinePipelineLayout, nullptr);
    }

    void VulkanRenderingSystem::DestroyDescriptors(VkDevice device){
        if (m_descriptorPool) vkDestroyDescriptorPool(device, m_descriptorPool, nullptr);
        if (m_transparentPool) vkDestroyDescriptorPool(device, m_transparentPool, nullptr);
        for (int i = 0; i < 2; i++) {
            if (m_lightPools[i])   vkDestroyDescriptorPool(device, m_lightPools[i], nullptr);
            if (m_svgfTemporalPools[i]) vkDestroyDescriptorPool(device, m_svgfTemporalPools[i], nullptr);
            if (m_svgfATrousPools[i]) vkDestroyDescriptorPool(device, m_svgfATrousPools[i], nullptr);
            if (m_taaPools[i])     vkDestroyDescriptorPool(device, m_taaPools[i], nullptr);
            if (m_blitPools[i])    vkDestroyDescriptorPool(device, m_blitPools[i], nullptr);
        }

        if (m_descriptorSetLayout) vkDestroyDescriptorSetLayout(device, m_descriptorSetLayout, nullptr);
        if (m_transparentSetLayout) vkDestroyDescriptorSetLayout(device, m_transparentSetLayout, nullptr);
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
        if (width == m_width && height == m_height) return;
        vkDeviceWaitIdle(m_device->GetNativeDevice());

        m_width = width;
        m_height = height;

        DestroyAttachments();
        CreateAttachments();

        UpdateDescriptorSets();

        InvalidateTemporalHistory();
    }

    void VulkanRenderingSystem::InvalidateTemporalHistory() {
        // Force temporal passes (TAA / SVGF) to re-bootstrap from current frame.
        m_frameIndex = 0;
        m_prevModelMats.clear();
        m_isFirstFrameRun = true;
        m_firstFrame = true;
    }

    void VulkanRenderingSystem::UpdateLightingAccelerationStructureDescriptor(uint32_t ping, VkAccelerationStructureKHR tlasHandle) {
        if (m_device == nullptr || tlasHandle == VK_NULL_HANDLE) {
            return;
        }

        VkWriteDescriptorSetAccelerationStructureKHR asWrite{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET_ACCELERATION_STRUCTURE_KHR};
        asWrite.accelerationStructureCount = 1;
        asWrite.pAccelerationStructures = &tlasHandle;

        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = m_lightingSets[ping];
        write.dstBinding = 4;
        write.dstArrayElement = 0;
        write.descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
        write.descriptorCount = 1;
        write.pNext = &asWrite;

        vkUpdateDescriptorSets(m_device->GetNativeDevice(), 1, &write, 0, nullptr);
    }

    void VulkanRenderingSystem::SetViewport(int x, int y, uint32_t width, uint32_t height) {
        m_viewX = x;
        m_viewY = y;
        m_viewW = width;
        m_viewH = height;
    }

    void VulkanRenderingSystem::SetEditorOverlayData(const EditorViewportOverlayData& overlayData) {
        m_editorOverlayData = overlayData;
        m_hasEditorOverlayData = true;
    }

    void VulkanRenderingSystem::UpdateBindlessTextureDescriptor(uint32_t textureIndex) {
        if (m_device == nullptr || m_descriptorSet == VK_NULL_HANDLE || textureIndex >= m_globalImageInfos.size()) {
            return;
        }

        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = m_descriptorSet;
        write.dstBinding = 0;
        write.dstArrayElement = textureIndex;
        write.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        write.descriptorCount = 1;
        write.pImageInfo = &m_globalImageInfos[textureIndex];
        vkUpdateDescriptorSets(m_device->GetNativeDevice(), 1, &write, 0, nullptr);
    }

    uint32_t VulkanRenderingSystem::GetOrLoadTextureIndex(const std::string& texturePath) {
        if (texturePath.empty()) {
            return static_cast<uint32_t>(-1);
        }

        auto cached = m_globalTexturePathCache.find(texturePath);
        if (cached != m_globalTexturePathCache.end()) {
            return cached->second;
        }

        if (m_globalTextures.size() >= 1024) {
            throw std::runtime_error("Bindless texture array is full.");
        }

        auto texture = std::make_unique<VulkanTexture>(m_device, m_resourceCommandPool.get(), texturePath);

        VkDescriptorImageInfo info{};
        info.imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        info.imageView = texture->GetImageView();
        info.sampler = texture->GetSampler();

        const uint32_t textureIndex = static_cast<uint32_t>(m_globalTextures.size());
        m_globalTextures.push_back(std::move(texture));
        m_globalImageInfos.push_back(info);
        m_globalTexturePathCache[texturePath] = textureIndex;
        UpdateBindlessTextureDescriptor(textureIndex);
        return textureIndex;
    }

    void VulkanRenderingSystem::ApplyMaterialOverride(MaterialData& targetMaterial, const VulkanMaterialSlotOverride& overrideData) {
        const uint32_t overrideMask = overrideData.materialInstance.overrideMask;
        if ((overrideMask & Resource::MaterialOverride_BaseColor) != 0u) {
            for (int channel = 0; channel < 4; ++channel) {
                targetMaterial.baseColorFactor[channel] = overrideData.materialInstance.factors.baseColorFactor[channel];
            }
        }
        if ((overrideMask & Resource::MaterialOverride_Metallic) != 0u) {
            targetMaterial.metallicFactor = overrideData.materialInstance.factors.metallicFactor;
        }
        if ((overrideMask & Resource::MaterialOverride_Roughness) != 0u) {
            targetMaterial.roughnessFactor = overrideData.materialInstance.factors.roughnessFactor;
        }
        if ((overrideMask & Resource::MaterialOverride_AlbedoTex) != 0u) {
            if (overrideData.materialInstance.textures.albedoTex >= 0) {
                targetMaterial.albedoTex = overrideData.materialInstance.textures.albedoTex;
            }
        }
        if ((overrideMask & Resource::MaterialOverride_NormalTex) != 0u) {
            if (overrideData.materialInstance.textures.normalTex >= 0) {
                targetMaterial.normalTex = overrideData.materialInstance.textures.normalTex;
            }
        }
        if ((overrideMask & Resource::MaterialOverride_OrmTex) != 0u) {
            if (overrideData.materialInstance.textures.ormTex >= 0) {
                targetMaterial.ormTex = overrideData.materialInstance.textures.ormTex;
            }
        }
        if ((overrideMask & Resource::MaterialOverride_EmissiveTex) != 0u) {
            if (overrideData.materialInstance.textures.emissiveTex >= 0) {
                targetMaterial.emissiveTex = overrideData.materialInstance.textures.emissiveTex;
            }
        }

        auto applyTexturePath = [this](const std::string& path, int* targetIndex) {
            if (targetIndex == nullptr) {
                return;
            }

            if (path.empty()) {
                return;
            }

            try {
                *targetIndex = static_cast<int>(GetOrLoadTextureIndex(path));
            } catch (const std::exception& exception) {
                std::cerr << "[VulkanRenderingSystem] WARNING: Failed to load override texture '" << path
                          << "'. Reason: " << exception.what() << std::endl;
            }
        };

        const VulkanTextureOverrideSet& textureOverrides = overrideData.textureOverrides;
        if (!textureOverrides.albedoTexturePath.empty()) {
            applyTexturePath(textureOverrides.albedoTexturePath, &targetMaterial.albedoTex);
        }
        if (!textureOverrides.normalTexturePath.empty()) {
            applyTexturePath(textureOverrides.normalTexturePath, &targetMaterial.normalTex);
        }
        if (!textureOverrides.ormTexturePath.empty()) {
            applyTexturePath(textureOverrides.ormTexturePath, &targetMaterial.ormTex);
        }
        if (!textureOverrides.emissiveTexturePath.empty()) {
            applyTexturePath(textureOverrides.emissiveTexturePath, &targetMaterial.emissiveTex);
        }
    }

    uint64_t VulkanRenderingSystem::ResolveMaterialBufferAddress(Entity entity, const VulkanModelComponent& modelComp, const VulkanModelResource& modelResource) {
        if (modelComp.GetMaterialOverrides().empty() || modelResource.materialsCpu.empty()) {
            return modelResource.matAddr;
        }

        const uint32_t entityId = static_cast<uint32_t>(entity);
        OverrideMaterialBufferCacheEntry& cacheEntry = m_overrideMaterialBuffers[entityId];
        if (cacheEntry.componentRevision == modelComp.GetResourceRevision()
            && cacheEntry.modelAssetPath == modelComp.getModelAssetPath()
            && cacheEntry.materialBufferAddress != 0) {
            return cacheEntry.materialBufferAddress;
        }

        std::vector<MaterialData> resolvedMaterials = modelResource.materialsCpu;
        for (const VulkanMaterialSlotOverride& overrideData : modelComp.GetMaterialOverrides()) {
            if (!overrideData.enabled || resolvedMaterials.empty()) {
                continue;
            }

            uint32_t targetMaterialIndex = overrideData.materialSlotIndex;
            if (overrideData.meshAssetIndex < modelResource.meshAssets.size()) {
                targetMaterialIndex = modelResource.meshAssets[overrideData.meshAssetIndex].materialIndex;
            }

            if (targetMaterialIndex >= resolvedMaterials.size()) {
                continue;
            }

            MaterialData targetMaterial = resolvedMaterials[targetMaterialIndex];
            if (overrideData.materialInstance.baseMaterialIndex < resolvedMaterials.size()
                && (overrideData.materialInstance.baseMaterialIndex != 0u || targetMaterialIndex == 0u)) {
                targetMaterial = resolvedMaterials[overrideData.materialInstance.baseMaterialIndex];
            }

            ApplyMaterialOverride(targetMaterial, overrideData);
            resolvedMaterials[targetMaterialIndex] = targetMaterial;
        }

        cacheEntry.materialBuffer = CreateBDABuffer(resolvedMaterials);
        cacheEntry.materialBufferAddress = cacheEntry.materialBuffer ? cacheEntry.materialBuffer->GetDeviceAddress() : modelResource.matAddr;
        cacheEntry.componentRevision = modelComp.GetResourceRevision();
        cacheEntry.modelAssetPath = modelComp.getModelAssetPath();
        return cacheEntry.materialBufferAddress;
    }

    VulkanModelResource& VulkanRenderingSystem::GetOrLoadModel(const std::string& path) {
        if (m_modelCache.find(path) != m_modelCache.end()) {
            return m_modelCache[path];
        }

        std::cout << "[RenderingSystem] Loading new model to GPU: " << path << std::endl;
        
        uint32_t currentTexOffset = static_cast<uint32_t>(m_globalTextures.size());
        MeshletModelBuilder builder;
        if (!builder.LoadAndSliceModel(path, currentTexOffset)) {
            throw std::runtime_error("Failed to load model asset: " + path);
        }

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
        res.materialCount = static_cast<uint32_t>(builder.GetMaterials().size());
        res.materialsCpu = builder.GetMaterials();
        res.materialAssets = builder.GetMaterialAssets();
        res.meshAssets = builder.GetMeshAssets();

        res.vertexBuffer   = CreateBDABuffer(builder.GetVertices());
        res.meshletBuffer  = CreateBDABuffer(builder.GetMeshlets());
        res.indexBuffer    = CreateBDABuffer(builder.GetMicroIndices());
        res.boundsBuffer   = CreateBDABuffer(builder.GetBounds());
        res.materialBuffer = CreateBDABuffer(builder.GetMaterials());

        res.vAddr   = res.vertexBuffer ? res.vertexBuffer->GetDeviceAddress() : 0;
        res.mAddr   = res.meshletBuffer ? res.meshletBuffer->GetDeviceAddress() : 0;
        res.iAddr   = res.indexBuffer ? res.indexBuffer->GetDeviceAddress() : 0;
        res.bAddr   = res.boundsBuffer ? res.boundsBuffer->GetDeviceAddress() : 0;
        res.matAddr = res.materialBuffer ? res.materialBuffer->GetDeviceAddress() : 0;

        const auto& vertices = builder.GetVertices();
        const auto& microIndices = builder.GetMicroIndices();
        const auto& meshlets = builder.GetMeshlets();

        std::vector<uint32_t> globalIndices;
        std::vector<uint32_t> primitiveMaterialIds;
        globalIndices.reserve(microIndices.size());
        primitiveMaterialIds.reserve(microIndices.size() / 3);
        for (const auto& m : meshlets) {
            for (uint32_t i = 0; i < m.triangleCount * 3; i++) {
                globalIndices.push_back(m.vertexOffset + microIndices[m.triangleOffset + i]); 
            }
            for (uint32_t tri = 0; tri < m.triangleCount; ++tri) {
                primitiveMaterialIds.push_back(m.materialID);
            }
        }

        uint32_t vertexCount = static_cast<uint32_t>(vertices.size());
        uint32_t indexCount = static_cast<uint32_t>(globalIndices.size());
        uint32_t vertexStride = vertices.empty() ? 0 : static_cast<uint32_t>(sizeof(vertices[0]));

        if (vertexCount > 0 && indexCount > 0) {
            res.globalIndexBuffer = CreateBDABuffer(globalIndices);
            res.globalIAddr = res.globalIndexBuffer->GetDeviceAddress();
            res.primitiveMaterialIdBuffer = CreateBDABuffer(primitiveMaterialIds);
            res.primMatIdAddr = res.primitiveMaterialIdBuffer ? res.primitiveMaterialIdBuffer->GetDeviceAddress() : 0;
            if (!res.primitiveMaterialIdBuffer) {
                res.materialCount = 0;
            }

            std::cout << "[VulkanBLAS] Triggering BLAS build for: " << path << std::endl;
            // BDA for globalIAddr
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
            write.dstArrayElement = currentTexOffset; 
            write.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            write.descriptorCount = res.textureCount;
            write.pImageInfo = m_globalImageInfos.data() + currentTexOffset;
            vkUpdateDescriptorSets(m_device->GetNativeDevice(), 1, &write, 0, nullptr);
        }

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
        struct PushConstants {
            uint64_t frameDataAddr;   
            uint64_t modelDataAddr; 
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
            uint64_t pointLightsAddr; 
            uint32_t stepSize;
            uint32_t pointLightCount;
        };

        struct LightingSpecializationData {
            int32_t giQualityLevel; 
            int32_t shadowQualityLevel;
        };

        VkPushConstantRange graphicsPushRange{}; 
        graphicsPushRange.stageFlags = VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT; 
        graphicsPushRange.size = sizeof(PushConstants); 

        VkPushConstantRange lightPushRange{}; 
        lightPushRange.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT; 
        lightPushRange.size = sizeof(LightingPushConstants);

        // global bindless texture
        VulkanDescriptorBuilder bindlessBuilder;
        VkDescriptorImageInfo dummyInfo = { m_gBufferSampler, m_gAlbedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };

        bindlessBuilder.BindImageArray(0, m_globalImageInfos.empty() ? &dummyInfo : m_globalImageInfos.data(), 1024, 
                                    m_globalImageInfos.empty() ? 1 : static_cast<uint32_t>(m_globalImageInfos.size()), 
                                    VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, VK_SHADER_STAGE_FRAGMENT_BIT, true);

        bindlessBuilder.Build(m_device, m_descriptorSetLayout, m_descriptorPool, m_descriptorSet);

        VkDescriptorSetLayoutBinding transparentBindings[1] = {};
        transparentBindings[0].binding = 0;
        transparentBindings[0].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        transparentBindings[0].descriptorCount = 1;
        transparentBindings[0].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;

        VkDescriptorSetLayoutCreateInfo transparentLayoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        transparentLayoutInfo.bindingCount = 1;
        transparentLayoutInfo.pBindings = transparentBindings;
        vkCreateDescriptorSetLayout(device, &transparentLayoutInfo, nullptr, &m_transparentSetLayout);

        VkDescriptorPoolSize transparentPoolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1};
        VkDescriptorPoolCreateInfo transparentPoolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
        transparentPoolInfo.poolSizeCount = 1;
        transparentPoolInfo.pPoolSizes = &transparentPoolSize;
        transparentPoolInfo.maxSets = 1;
        vkCreateDescriptorPool(device, &transparentPoolInfo, nullptr, &m_transparentPool);

        VkDescriptorSetAllocateInfo transparentAllocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        transparentAllocInfo.descriptorPool = m_transparentPool;
        transparentAllocInfo.descriptorSetCount = 1;
        transparentAllocInfo.pSetLayouts = &m_transparentSetLayout;
        vkAllocateDescriptorSets(device, &transparentAllocInfo, &m_transparentSet);

        for (uint32_t ping = 0; ping < 2; ping++) {
            VkDevice dev = m_device->GetNativeDevice();

            // A. Lighting Pass 
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

            // B. SVGF Temporal Pass
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

            // C. SVGF A-Trous Pass
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

            // D. TAA Pass
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

            // E. Blit Pass
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

        // Mesh Shaders
        VkShaderModule taskShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "car_task.spv"));
        VkShaderModule meshShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "car_mesh.spv"));
        VkShaderModule fragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "car_frag.spv"));
        VkShaderModule transparentFragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "transparent_frag.spv"));

        // Fullscreen Quad 
        VkShaderModule fsVertShader    = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "lighting_vert.spv"));
        VkShaderModule lightFragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "lighting_uber.spv"));
        VkShaderModule svgfTemporalFragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "svgf_temporal_frag.spv"));
        VkShaderModule svgfATrousFragShader   = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "svgf_atrous_frag.spv"));
        VkShaderModule taaFragShader   = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "taa_frag.spv"));
        VkShaderModule blitFragShader  = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "blit_frag.spv"));
   
        // Graphics Pipeline (G-Buffer)
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

        VkDescriptorSetLayout transparentSetLayouts[2] = { m_transparentSetLayout, m_descriptorSetLayout };
        VkPipelineLayoutCreateInfo transparentPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        transparentPipeLayoutInfo.setLayoutCount = 2;
        transparentPipeLayoutInfo.pSetLayouts = transparentSetLayouts;
        transparentPipeLayoutInfo.pushConstantRangeCount = 1;
        transparentPipeLayoutInfo.pPushConstantRanges = &graphicsPushRange;
        vkCreatePipelineLayout(device, &transparentPipeLayoutInfo, nullptr, &m_transparentPipelineLayout);

        m_transparentPipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_TASK_BIT_EXT, taskShader)
            .AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, transparentFragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(true)
            .SetPipelineLayout(m_transparentPipelineLayout)
            .Build(m_device, { m_renderer->GetSwapchainFormat() }, VK_FORMAT_D32_SFLOAT);


        LightingSpecializationData specData{};
        specData.giQualityLevel = m_lightingProfile.giQualityLevel;
        specData.shadowQualityLevel = m_lightingProfile.shadowQualityLevel;

        VkSpecializationMapEntry specEntries[2] = {};
        // layout(constant_id = 0) const int GI_QUALITY_LEVEL;
        specEntries[0].constantID = 0;
        specEntries[0].offset = offsetof(LightingSpecializationData, giQualityLevel);
        specEntries[0].size = sizeof(int32_t);

        // layout(constant_id = 1) const int ENABLE_SHADOWS;
        specEntries[1].constantID = 1;
        specEntries[1].offset = offsetof(LightingSpecializationData, shadowQualityLevel);
        specEntries[1].size = sizeof(int32_t);

        VkSpecializationInfo specInfo{};
        specInfo.mapEntryCount = 2;
        specInfo.pMapEntries = specEntries;
        specInfo.dataSize = sizeof(LightingSpecializationData);
        specInfo.pData = &specData;

        // Lighting Pipeline
        VkDescriptorSetLayout lightSetLayouts[2] = { m_lightingLayouts[0], m_descriptorSetLayout };
        VkPipelineLayoutCreateInfo lightPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        lightPipeLayoutInfo.setLayoutCount = 2; 
        lightPipeLayoutInfo.pSetLayouts = lightSetLayouts; 
        lightPipeLayoutInfo.pushConstantRangeCount = 1; 
        lightPipeLayoutInfo.pPushConstantRanges = &lightPushRange; 
        vkCreatePipelineLayout(device, &lightPipeLayoutInfo, nullptr, &m_lightingPipelineLayout);

        m_lightingPipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, fsVertShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, lightFragShader, &specInfo)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false).AddColorBlendAttachment(false) // DirectLight, NoisyGI
            .SetPipelineLayout(m_lightingPipelineLayout)
            .Build(m_device, { VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED);

        // SVGF Temporal Pipeline layout GI History & Moments
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
            .AddColorBlendAttachment(false) 
            .AddColorBlendAttachment(false) 
            .SetPipelineLayout(m_svgfTemporalPipelineLayout)
            .Build(m_device, { VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16_SFLOAT }, VK_FORMAT_UNDEFINED);

        // SVGF A-Trous Pipeline 
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
            .AddColorBlendAttachment(false) // only DenoisedGI
            .SetPipelineLayout(m_svgfATrousPipelineLayout)
            .Build(m_device, { VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED);

        // 4.5 TAA Pipeline
        VkPipelineLayoutCreateInfo taaPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        taaPipeLayoutInfo.setLayoutCount = 1; 
        taaPipeLayoutInfo.pSetLayouts = &m_taaLayouts[0]; 
        taaPipeLayoutInfo.pushConstantRangeCount = 1;             
        taaPipeLayoutInfo.pPushConstantRanges = &lightPushRange;
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

        // Blit Pipeline 
        VkPipelineLayoutCreateInfo blitPipeLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; 
        blitPipeLayoutInfo.setLayoutCount = 1; 
        blitPipeLayoutInfo.pSetLayouts = &m_blitLayouts[0]; 
        blitPipeLayoutInfo.pushConstantRangeCount = 1;             
        blitPipeLayoutInfo.pPushConstantRanges = &lightPushRange;
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
            .Build(m_device, { m_renderer->GetSwapchainFormat() }, VK_FORMAT_D32_SFLOAT);

        vkDestroyShaderModule(device, taskShader, nullptr);
        vkDestroyShaderModule(device, meshShader, nullptr);
        vkDestroyShaderModule(device, fragShader, nullptr);
        vkDestroyShaderModule(device, transparentFragShader, nullptr);
        vkDestroyShaderModule(device, fsVertShader, nullptr);
        vkDestroyShaderModule(device, lightFragShader, nullptr);
        vkDestroyShaderModule(device, svgfTemporalFragShader, nullptr);
        vkDestroyShaderModule(device, svgfATrousFragShader, nullptr);
        vkDestroyShaderModule(device, taaFragShader, nullptr);
        vkDestroyShaderModule(device, blitFragShader, nullptr);

        CreateEditorGridPipeline();
        CreateDebugLinePipeline();
        std::cout << "[RenderingSystem] Pipelines Built Successfully." << std::endl;
    }

    void VulkanRenderingSystem::CreateEditorGridPipeline() {
        VkDevice device = m_device->GetNativeDevice();
        const std::string SHADER_DIR = "C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/shader/";

        VkShaderModule vertShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "editor_grid_vert.spv"));
        VkShaderModule fragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "editor_grid_frag.spv"));

        VkPipelineShaderStageCreateInfo shaderStages[] = {
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_VERTEX_BIT, vertShader, "main", nullptr},
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_FRAGMENT_BIT, fragShader, "main", nullptr}
        };

        VkPipelineVertexInputStateCreateInfo vertexInputInfo{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};

        VkPipelineInputAssemblyStateCreateInfo inputAssembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
        inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
        inputAssembly.primitiveRestartEnable = VK_FALSE;

        VkPipelineViewportStateCreateInfo viewportState{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
        viewportState.viewportCount = 1;
        viewportState.scissorCount = 1;

        VkPipelineRasterizationStateCreateInfo rasterizer{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
        rasterizer.lineWidth = 1.0f;
        rasterizer.cullMode = VK_CULL_MODE_NONE;
        rasterizer.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;

        VkPipelineMultisampleStateCreateInfo multisampling{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
        multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

        VkPipelineDepthStencilStateCreateInfo depthStencil{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
        depthStencil.depthTestEnable = VK_FALSE;
        depthStencil.depthWriteEnable = VK_FALSE;

        VkPipelineColorBlendAttachmentState colorBlendAttachment{};
        colorBlendAttachment.colorWriteMask = 0xF;
        colorBlendAttachment.blendEnable = VK_TRUE;
        colorBlendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
        colorBlendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        colorBlendAttachment.colorBlendOp = VK_BLEND_OP_ADD;
        colorBlendAttachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
        colorBlendAttachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        colorBlendAttachment.alphaBlendOp = VK_BLEND_OP_ADD;

        VkPipelineColorBlendStateCreateInfo colorBlending{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
        colorBlending.attachmentCount = 1;
        colorBlending.pAttachments = &colorBlendAttachment;

        VkDynamicState dynamicStates[] = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };
        VkPipelineDynamicStateCreateInfo dynamicState{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
        dynamicState.dynamicStateCount = 2;
        dynamicState.pDynamicStates = dynamicStates;

        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
        pushRange.offset = 0;
        pushRange.size = sizeof(EditorGridPushConstants);

        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(device, &pipelineLayoutInfo, nullptr, &m_editorGridPipelineLayout);

        VkFormat colorFormat = m_renderer->GetSwapchainFormat();
        VkPipelineRenderingCreateInfo renderingInfo{VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO};
        renderingInfo.colorAttachmentCount = 1;
        renderingInfo.pColorAttachmentFormats = &colorFormat;
        renderingInfo.depthAttachmentFormat = VK_FORMAT_D32_SFLOAT;

        VkGraphicsPipelineCreateInfo pipelineInfo{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        pipelineInfo.pNext = &renderingInfo;
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = shaderStages;
        pipelineInfo.pVertexInputState = &vertexInputInfo;
        pipelineInfo.pInputAssemblyState = &inputAssembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &rasterizer;
        pipelineInfo.pMultisampleState = &multisampling;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &colorBlending;
        pipelineInfo.pDynamicState = &dynamicState;
        pipelineInfo.layout = m_editorGridPipelineLayout;

        if (vkCreateGraphicsPipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &m_editorGridPipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create editor grid pipeline!");
        }

        vkDestroyShaderModule(device, vertShader, nullptr);
        vkDestroyShaderModule(device, fragShader, nullptr);
    }

    void VulkanRenderingSystem::CreateDebugLinePipeline() {
        VkDevice device = m_device->GetNativeDevice();
        const std::string SHADER_DIR = "C:/Lizeral Engine/LizeralEngine0.0.1/engine/source/shader/";

        VkShaderModule vertShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "debug_line_vert.spv"));
        VkShaderModule fragShader = CreateShaderModule(device, ReadShaderFile(SHADER_DIR + "debug_line_frag.spv"));

        VkPipelineShaderStageCreateInfo shaderStages[] = {
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_VERTEX_BIT, vertShader, "main", nullptr},
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_FRAGMENT_BIT, fragShader, "main", nullptr}
        };

        // Vertex Input 
        VkVertexInputBindingDescription bindingDesc{};
        bindingDesc.binding = 0;
        bindingDesc.stride = sizeof(DebugLineVertex);
        bindingDesc.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;

        VkVertexInputAttributeDescription attrDesc[2] = {};
        attrDesc[0].binding = 0; attrDesc[0].location = 0; attrDesc[0].format = VK_FORMAT_R32G32B32_SFLOAT; attrDesc[0].offset = offsetof(DebugLineVertex, position);
        attrDesc[1].binding = 0; attrDesc[1].location = 1; attrDesc[1].format = VK_FORMAT_R32G32B32_SFLOAT; attrDesc[1].offset = offsetof(DebugLineVertex, color);

        VkPipelineVertexInputStateCreateInfo vertexInputInfo{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
        vertexInputInfo.vertexBindingDescriptionCount = 1; vertexInputInfo.pVertexBindingDescriptions = &bindingDesc;
        vertexInputInfo.vertexAttributeDescriptionCount = 2; vertexInputInfo.pVertexAttributeDescriptions = attrDesc;

        // Topology = LINE_LIST
        VkPipelineInputAssemblyStateCreateInfo inputAssembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
        inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_LINE_LIST;
        inputAssembly.primitiveRestartEnable = VK_FALSE;

        VkPipelineViewportStateCreateInfo viewportState{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
        viewportState.viewportCount = 1; viewportState.scissorCount = 1;

        VkPipelineRasterizationStateCreateInfo rasterizer{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
        rasterizer.lineWidth = 1.0f; 
        rasterizer.cullMode = VK_CULL_MODE_NONE;
        rasterizer.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;

        VkPipelineMultisampleStateCreateInfo multisampling{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
        multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

        VkPipelineDepthStencilStateCreateInfo depthStencil{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
        depthStencil.depthTestEnable = VK_FALSE; depthStencil.depthWriteEnable = VK_FALSE;

        VkPipelineColorBlendAttachmentState colorBlendAttachment{};
        colorBlendAttachment.colorWriteMask = 0xF; colorBlendAttachment.blendEnable = VK_FALSE;
        VkPipelineColorBlendStateCreateInfo colorBlending{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
        colorBlending.attachmentCount = 1; colorBlending.pAttachments = &colorBlendAttachment;

        VkDynamicState dynamicStates[] = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };
        VkPipelineDynamicStateCreateInfo dynamicState{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
        dynamicState.dynamicStateCount = 2; dynamicState.pDynamicStates = dynamicStates;

        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_VERTEX_BIT; pushRange.offset = 0; pushRange.size = sizeof(Matrix4x4);

        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.pushConstantRangeCount = 1; pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(device, &pipelineLayoutInfo, nullptr, &m_debugLinePipelineLayout);

        VkFormat colorFormat = m_renderer->GetSwapchainFormat();
        VkPipelineRenderingCreateInfo renderingInfo{VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO};
        renderingInfo.colorAttachmentCount = 1; renderingInfo.pColorAttachmentFormats = &colorFormat;
        renderingInfo.depthAttachmentFormat = VK_FORMAT_D32_SFLOAT;

        VkGraphicsPipelineCreateInfo pipelineInfo{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        pipelineInfo.pNext = &renderingInfo;
        pipelineInfo.stageCount = 2; pipelineInfo.pStages = shaderStages;
        pipelineInfo.pVertexInputState = &vertexInputInfo; pipelineInfo.pInputAssemblyState = &inputAssembly;
        pipelineInfo.pViewportState = &viewportState; pipelineInfo.pRasterizationState = &rasterizer;
        pipelineInfo.pMultisampleState = &multisampling; pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &colorBlending; pipelineInfo.pDynamicState = &dynamicState;
        pipelineInfo.layout = m_debugLinePipelineLayout;

        if (vkCreateGraphicsPipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &m_debugLinePipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create debug line pipeline!");
        }

        vkDestroyShaderModule(device, vertShader, nullptr);
        vkDestroyShaderModule(device, fragShader, nullptr);
    }


    void VulkanRenderingSystem::Tick(Registry& registry, float deltaTime,const std::vector<DebugLineVertex>& debugLines) {
        if (!m_CmdDrawMeshTasksEXT) {
            m_CmdDrawMeshTasksEXT = (PFN_vkCmdDrawMeshTasksEXT)vkGetDeviceProcAddr(m_device->GetNativeDevice(), "vkCmdDrawMeshTasksEXT");
        }

        const EditorViewportOverlayData* activeOverlay = m_hasEditorOverlayData ? &m_editorOverlayData : nullptr;
        const std::vector<DebugLineVertex>& activeDebugLines =
            activeOverlay ? activeOverlay->worldLines : debugLines;

        // Keep render targets in sync with the real swapchain size.
        // Some fullscreen / HiDPI transitions may bypass widget resize callbacks.
        VkExtent2D swapExtent = m_renderer->GetSwapchainExtent();
        if (swapExtent.width > 0 && swapExtent.height > 0 &&
            (swapExtent.width != m_width || swapExtent.height != m_height)) {
            Resize(swapExtent.width, swapExtent.height);
        }

        VkCommandBuffer cmd = m_renderer->BeginFrame();
        if (cmd == VK_NULL_HANDLE) return;

        // Advance temporal frame index only when a frame is actually rendered.
        m_frameIndex++;

        // TAA Jitter caculate
        uint32_t jitterIndex = m_frameIndex % 8 + 1; 
        float jitterX = CreateHaltonSequence(jitterIndex, 2) - 0.5f; 
        float jitterY = CreateHaltonSequence(jitterIndex, 3) - 0.5f;
        float ndcJitterX = (jitterX * 2.0f) / static_cast<float>(m_width > 0 ? m_width : 1); 
        float ndcJitterY = (jitterY * 2.0f) / static_cast<float>(m_height > 0 ? m_height : 1);

        // Descriptor/history ping-pong must follow renderer frame slot (fence-protected),
        // not temporal frame index (which can be reset by editor interactions).
        uint32_t ping = m_renderer->GetCurrentFrameSlot();

        VkViewport viewport{ 0.0f, 0.0f, (float)m_width, (float)m_height, 0.0f, 1.0f };
        VkRect2D scissor{ {0, 0}, {m_width, m_height} };

        Matrix4x4 viewMat, projMatUnjittered, projMatJittered;
        Vector3 cameraPos;
        auto cameraView = registry.view<TransformComponent, CameraComponent>();
        
        for (auto entity : cameraView) {
            auto& camera = cameraView.get<CameraComponent>(entity);
            auto& transform = cameraView.get<TransformComponent>(entity);
            cameraPos = transform.getPosition();

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

        // No implicit directional sun by default.
        // Directional lighting is enabled only when a DirectionLightComponent exists.
        Lizeral::Vector3 lightDir(0.0f, -1.0f, 0.0f);
        Lizeral::Vector3 lightColor(1.0f, 1.0f, 1.0f);
        float lightIntensity = 0.0f;

        auto lightView = registry.view<TransformComponent, DirectionLightComponent>();
        for (auto entity : lightView) {
            auto& trans = lightView.get<TransformComponent>(entity);
            auto& light = lightView.get<DirectionLightComponent>(entity);
            lightDir = trans.getForward().normalize();
            lightColor = light.getColor(); 
            lightIntensity = light.getIntensity();
            break; 
        }

        Matrix4x4 currentVP = projMatUnjittered * viewMat; //projMatUnjittered is important !!!
        if (m_isFirstFrameRun) {
            m_prevViewProj = currentVP;
        }

        GlobalFrameData frameData{};
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

        RTInstanceDesc* mappedDesc = static_cast<RTInstanceDesc*>(m_rtInstanceBuffer->Map());

        std::vector<VkAccelerationStructureInstanceKHR> tlasInstances;
        uint32_t customInstanceId = 0;

        auto modelView = registry.view<TransformComponent, VulkanModelComponent>();
        struct TransparentDrawItem {
            Entity entity = null_entity;
            float distanceSq = 0.0f;
        };
        std::vector<TransparentDrawItem> transparentDraws;
        std::unordered_map<uint32_t, uint32_t> entityToInstanceIndex;

        auto HasTransparentMaterial = [](const VulkanModelResource& resource) {
            for (const auto& material : resource.materialsCpu) {
                if (material.alphaMode == static_cast<int>(MaterialAlphaMode::Blend)) {
                    return true;
                }
            }
            return false;
        };

        for (auto entity : modelView) {
            auto& transform = modelView.get<TransformComponent>(entity);
            auto& modelComp = modelView.get<VulkanModelComponent>(entity);
            
            if (!modelComp.IsVisible() || !modelComp.HasModelAsset()) continue;

            auto& res = GetOrLoadModel(modelComp.getModelAssetPath());
            if (!res.IsValid() || !res.blas) continue;
            const uint64_t resolvedMaterialBuffer = ResolveMaterialBufferAddress(entity, modelComp, res);

            mappedDesc[customInstanceId].vertexBuffer = res.vAddr;
            mappedDesc[customInstanceId].indexBuffer  = res.globalIAddr;
            mappedDesc[customInstanceId].materialBuffer = resolvedMaterialBuffer;
            mappedDesc[customInstanceId].primitiveMaterialIdBuffer = res.primMatIdAddr;
            mappedDesc[customInstanceId].textureOffset  = res.textureOffset;
            mappedDesc[customInstanceId].materialCount = res.materialCount;
            mappedDesc[customInstanceId].padding[0] = 0;
            mappedDesc[customInstanceId].padding[1] = 0;

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

            if (HasTransparentMaterial(res)) {
                TransparentDrawItem item;
                item.entity = entity;
                item.distanceSq = transform.getPosition().squaredDistance(cameraPos);
                transparentDraws.push_back(item);
            }
        }
        m_rtInstanceBuffer->Unmap();

        std::sort(transparentDraws.begin(), transparentDraws.end(), [](const TransparentDrawItem& lhs, const TransparentDrawItem& rhs) {
            return lhs.distanceSq > rhs.distanceSq;
        });

        if (!tlasInstances.empty()) {
            m_tlas->Build(cmd, ping, tlasInstances, false); 
            VkMemoryBarrier memoryBarrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
            memoryBarrier.srcAccessMask = VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR;
            memoryBarrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
            vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
        }

        // hot swap Descriptor
        VkAccelerationStructureKHR currentTlas = m_tlas->GetHandle(ping);
        UpdateLightingAccelerationStructureDescriptor(ping, currentTlas);

        //  G-Buffer Pass
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

        struct PushConstants {
            uint64_t frameDataAddr; 
            uint64_t modelDataAddr;   
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
            if (!modelComp.IsVisible() || !modelComp.HasModelAsset()) continue;

            const auto& res = GetOrLoadModel(modelComp.getModelAssetPath());
            if (!res.IsValid()) continue;
            const uint64_t resolvedMaterialBuffer = ResolveMaterialBufferAddress(entity, modelComp, res);

            Matrix4x4 currentModel = transform.getMatrix();
            uint32_t entityId = static_cast<uint32_t>(entity);
            Matrix4x4 prevModel = m_isFirstFrameRun ? currentModel : m_prevModelMats[entityId];

            PushConstants pushData{};

            if (currentInstanceIndex >= 10000) {
                std::cerr << "WARNING: Max instances (10000) reached!" << std::endl;
                break;
            }

            mappedInstanceData[currentInstanceIndex].Model = currentModel.transpose();
            mappedInstanceData[currentInstanceIndex].prevModel = prevModel.transpose();
            entityToInstanceIndex[static_cast<uint32_t>(entity)] = currentInstanceIndex;
            
            pushData.frameDataAddr = m_frameResource.addr;
            pushData.modelDataAddr = baseInstanceAddr + (currentInstanceIndex * sizeof(GPUInstanceData));
            pushData.vertexBuffer = res.vAddr; pushData.meshletBuffer = res.mAddr; pushData.indexBuffer = res.iAddr; 
            pushData.boundsBuffer = res.bAddr; pushData.materialBuffer = resolvedMaterialBuffer;
            pushData.totalMeshlets = res.totalMeshlets; pushData.textureOffset = 0; 

            vkCmdPushConstants(cmd, m_graphicsPipelineLayout, VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(PushConstants), &pushData);
            m_prevModelMats[entityId] = currentModel;
            m_CmdDrawMeshTasksEXT(cmd, (res.totalMeshlets + 63) / 64, 1, 1);
            currentInstanceIndex++;
        }

        std::vector<GPUPointLight> activePointLights;
        auto plView = registry.view<TransformComponent, PointLightComponent>();
        for (auto entity : plView) {
            auto& trans = plView.get<TransformComponent>(entity);
            auto& pl = plView.get<PointLightComponent>(entity);
            
            activePointLights.push_back({
                Vector4(trans.getPosition().x, trans.getPosition().y, trans.getPosition().z, pl.getSpotRadius()),
                Vector4(pl.getSpotLightColor().x, pl.getSpotLightColor().y, pl.getSpotLightColor().z, pl.getSpotIntensity())
            });
        }


        uint32_t lightCount = static_cast<uint32_t>(activePointLights.size());
        uint32_t maxLightCapacity = 1024; 
        if (lightCount > maxLightCapacity) {
            std::cerr << "[Warning] Point lights count exceeds capacity (1024), truncating!" << std::endl;
            lightCount = maxLightCapacity;
        }
        if (lightCount > 0) { 
            m_pointLightBuffer->WriteData(activePointLights.data(), lightCount * sizeof(GPUPointLight));
        }

        vkCmdEndRendering(cmd);

        m_prevViewProj = currentViewProjUnjittered;
        if (m_isFirstFrameRun) {
            TransitionImageLayout(cmd, m_gDirectLight.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_gNoisyGI.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_gDenoisedGI.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_gDenoisedGITemp.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
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
            uint64_t pointLightsAddr;
            uint32_t stepSize;
            uint32_t pointLightCount;
        };

        LightingPushConstants lightPc{};
        lightPc.frameDataAddr = m_frameResource.addr;
        lightPc.instanceDescAddr = m_rtInstanceBuffer->GetDeviceAddress();
        lightPc.pointLightsAddr = m_pointLightBuffer ? m_pointLightBuffer->GetDeviceAddress() : 0;
        lightPc.pointLightCount = lightCount;


        auto DrawFullscreenPass = [&](VkPipeline pipeline, VkPipelineLayout layout, VkDescriptorSet set, GBufferAttachment* outputs, uint32_t outputCount, bool bindExtraSet = false) {
            VkRenderingAttachmentInfo attInfos[2] = {};
            for(uint32_t i=0; i<outputCount; i++) {
                TransitionImageLayout(cmd, outputs[i].image, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
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

        if (activeOverlay && activeOverlay->enabled && activeOverlay->grid.enabled && m_editorGridPipeline != VK_NULL_HANDLE) {
            const float minorSpacing = std::max(activeOverlay->grid.minorSpacing, 0.001f);
            const float majorSpacing = std::max(activeOverlay->grid.majorSpacing, minorSpacing);

            EditorGridPushConstants gridPc{};
            gridPc.invViewProj = currentVP.inverse().transpose();
            gridPc.cameraPosAndPlaneHeight = Vector4(cameraPos, activeOverlay->grid.planeHeight);
            gridPc.viewportSizeAndSpacing =
                Vector4(static_cast<float>(m_width), static_cast<float>(m_height), minorSpacing, majorSpacing);
            gridPc.fadeAndOpacity =
                Vector4(std::max(activeOverlay->grid.fadeDistance, 1.0f),
                        activeOverlay->grid.majorOpacity,
                        activeOverlay->grid.minorOpacity,
                        activeOverlay->grid.axisOpacity);

            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_editorGridPipeline);
            vkCmdPushConstants(cmd, m_editorGridPipelineLayout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(EditorGridPushConstants), &gridPc);
            vkCmdDraw(cmd, 3, 1, 0, 0);
        }

        if (!transparentDraws.empty()) {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_transparentPipeline);
            VkDescriptorSet transparentSets[2] = { m_transparentSet, m_descriptorSet };
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_transparentPipelineLayout, 0, 2, transparentSets, 0, nullptr);

            for (const TransparentDrawItem& draw : transparentDraws) {
                auto instanceIt = entityToInstanceIndex.find(static_cast<uint32_t>(draw.entity));
                if (instanceIt == entityToInstanceIndex.end()) {
                    continue;
                }

                auto& transform = modelView.get<TransformComponent>(draw.entity);
                auto& modelComp = modelView.get<VulkanModelComponent>(draw.entity);
                if (!modelComp.IsVisible() || !modelComp.HasModelAsset()) {
                    continue;
                }

                const auto& res = GetOrLoadModel(modelComp.getModelAssetPath());
                if (!res.IsValid()) {
                    continue;
                }

                PushConstants pushData{};
                pushData.frameDataAddr = m_frameResource.addr;
                pushData.modelDataAddr = baseInstanceAddr + (instanceIt->second * sizeof(GPUInstanceData));
                pushData.vertexBuffer = res.vAddr;
                pushData.meshletBuffer = res.mAddr;
                pushData.indexBuffer = res.iAddr;
                pushData.boundsBuffer = res.bAddr;
                pushData.materialBuffer = ResolveMaterialBufferAddress(draw.entity, modelComp, res);
                pushData.totalMeshlets = res.totalMeshlets;
                pushData.textureOffset = 0;

                vkCmdPushConstants(cmd, m_transparentPipelineLayout, VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(PushConstants), &pushData);
                m_CmdDrawMeshTasksEXT(cmd, (res.totalMeshlets + 63) / 64, 1, 1);
            }
        }

        // 6. Draw Lines
        if (!activeDebugLines.empty()) {
            size_t bufferSize = activeDebugLines.size() * sizeof(Lizeral::DebugLineVertex);
            if (bufferSize > m_maxDebugLineBufferSize) bufferSize = m_maxDebugLineBufferSize;

            m_debugLineBuffer->WriteData((void*)activeDebugLines.data(), bufferSize);

            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_debugLinePipeline);
            Matrix4x4 lineVP = currentVP.transpose(); 
            vkCmdPushConstants(cmd, m_debugLinePipelineLayout, VK_SHADER_STAGE_VERTEX_BIT, 0, sizeof(Matrix4x4), &lineVP);

            VkDeviceSize offsets[] = {0};
            VkBuffer rawBuffer = m_debugLineBuffer->GetNativeBuffer();
            vkCmdBindVertexBuffers(cmd, 0, 1, &rawBuffer, offsets);
            
            uint32_t drawCount = static_cast<uint32_t>(bufferSize / sizeof(Lizeral::DebugLineVertex));
            vkCmdDraw(cmd, drawCount, 1, 0, 0);
        }

        m_renderer->EndRendering(cmd);

        m_renderer->EndFrame();
        m_isFirstFrameRun = m_firstFrame = false;
    }


    void VulkanRenderingSystem::UpdateDescriptorSets() {
        VkDevice device = m_device->GetNativeDevice();

        VkDescriptorImageInfo transparentDepthInfo {
            m_gBufferSampler,
            m_gDepth.view,
            VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL
        };
        VkWriteDescriptorSet transparentWrite{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        transparentWrite.dstSet = m_transparentSet;
        transparentWrite.dstBinding = 0;
        transparentWrite.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        transparentWrite.descriptorCount = 1;
        transparentWrite.pImageInfo = &transparentDepthInfo;
        vkUpdateDescriptorSets(device, 1, &transparentWrite, 0, nullptr);

        VkDescriptorImageInfo gInfos[4] = {
            { m_gBufferSampler, m_gAlbedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { m_gBufferSampler, m_gNormalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { m_gBufferSampler, m_gDepth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { m_gBufferSampler, m_gVelocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
        };

        for (uint32_t i = 0; i < 2; i++) {
            uint32_t ping = i;
            uint32_t pong = (i + 1) % 2;

            // Lighting Pass Sets 
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
            if (m_tlas) {
                UpdateLightingAccelerationStructureDescriptor(ping, m_tlas->GetHandle(ping));
            }

            // SVGF Temporal Pass Sets
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

            // SVGF A-Trous Pass Sets 
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

            // TAA Pass Sets
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

            // Blit Pass Sets
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
