#pragma once

#include "editor/overlay/EditorViewportOverlayTypes.h"
#include "runtime/function/ecs/registry.h"
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"
#include "runtime/function/render/rhi/vulkan/VulkanBuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandPool.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandBuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanTexture.h"
#include "runtime/function/render/rhi/vulkan/VulkanTLAS.h"
#include "runtime/function/render/rhi/vulkan/VulkanBLAS.h"
#include "runtime/function/Vulkan_res_type/VulkanModelResource.h"
#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/core/math/matrix4.h"
#include "runtime/core/math/vector2.h"

#include <memory>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstdint>

namespace Lizeral {

    class VulkanModelComponent;
    struct VulkanMaterialSlotOverride;

    struct GBufferAttachment {
        VkImage image = VK_NULL_HANDLE;
        VmaAllocation allocation = VK_NULL_HANDLE;
        VkImageView view = VK_NULL_HANDLE;
        VkFormat format;
    };

    struct GlobalFrameData {
        Matrix4x4 viewProj;
        Matrix4x4 invViewProj;
        Matrix4x4 prevViewProj;
        Vector3 cameraPos;
        float padding1;
        Vector3 lightDir;
        float lightIntensity;
        Vector3 lightColor;
        uint32_t frameIndex;
        Vector2 jitter;
        Vector2 padding2; 
    };

    struct GPUPointLight {
        Vector4 posAndRadius;      // x, y, z, w(radius)
        Vector4 colorAndIntensity; // r, g, b, a(intensity)
    };

    // activePointLights.push_back({trans.getPosition(), pl.getRadius(), pl.getLightColor(), pl.getIntensity()});

    struct VulkanFrameResource {
        std::unique_ptr<VulkanBuffer> buffer;
        uint64_t addr = 0;
        bool IsValid() const { return addr != 0; }
    };

    struct RTInstanceDesc {
        uint64_t vertexBuffer;
        uint64_t indexBuffer;
        uint64_t materialBuffer;
        uint64_t primitiveMaterialIdBuffer;
        uint32_t textureOffset;
        uint32_t materialCount;
        uint32_t padding[2];
    };

    struct GPUInstanceData {
        Matrix4x4 Model;
        Matrix4x4 prevModel;
    };

    class VulkanRenderingSystem {
    public:
        struct LightingProfile {
            int32_t giQualityLevel = 2;      // 0: off, 1: SSGI, 2: RTGI
            int32_t shadowQualityLevel = 1;  // 0: hard, 1: soft
        };

        enum class RenderPipelinePreset : uint8_t {
            Stable = 0,   // GI off, shadow hard
            SSGI = 1,     // GI SSGI, shadow hard
            RTGI = 2      // GI RTGI, shadow soft
        };

        VulkanRenderingSystem() = default;
        ~VulkanRenderingSystem() {}

        void Initialize(VulkanContext* context, VulkanDevice* device, VulkanRenderer* renderer, uint32_t width, uint32_t height);
        
        void Tick(Registry& registry, float deltaTime, const std::vector<Lizeral::DebugLineVertex>& debugLines = {});
        
        void Shutdown();
        void Resize(uint32_t width, uint32_t height);
        void InvalidateTemporalHistory();
        void SetViewport(int x, int y, uint32_t width, uint32_t height);
        void SetEditorOverlayData(const EditorViewportOverlayData& overlayData);
        void SetRenderPipelinePreset(RenderPipelinePreset preset);
        RenderPipelinePreset GetRenderPipelinePreset() const { return m_renderPreset; }
        void SetLightingProfile(const LightingProfile& profile);
        LightingProfile GetLightingProfile() const { return m_lightingProfile; }
        void ResetLightingProfileToPreset();
        void SetDefaultFBO(unsigned int fbo) { /* VulkanRenderer handles target */ }

        void WaitIdle() {
            if (m_device) m_device->WaitIdle();
        }

    private:
        // RHI handle
        VulkanContext* m_context { nullptr };
        VulkanDevice* m_device { nullptr };
        VulkanRenderer* m_renderer { nullptr };

        std::unique_ptr<VulkanCommandPool> m_resourceCommandPool;
        std::unique_ptr<VulkanTLAS> m_tlas;

        // G-Buffer 
        GBufferAttachment m_gAlbedoMetallic;
        GBufferAttachment m_gNormalRoughness;
        GBufferAttachment m_gDepth;
        GBufferAttachment m_gVelocity;
        GBufferAttachment m_gDirectLight;
        GBufferAttachment m_gNoisyGI;
        GBufferAttachment m_gDenoisedGI;
        GBufferAttachment m_gDenoisedGITemp;
        GBufferAttachment m_gGIHistory[2];      
        GBufferAttachment m_gMomentsHistory[2];   
        GBufferAttachment m_gHistory[2];
        VkSampler m_gBufferSampler { VK_NULL_HANDLE };

        // VkDescriptors
        VkDescriptorSetLayout m_descriptorSetLayout { VK_NULL_HANDLE };
        VkDescriptorPool m_descriptorPool { VK_NULL_HANDLE };
        VkDescriptorSet m_descriptorSet { VK_NULL_HANDLE };

        VkDescriptorSetLayout m_lightingLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_lightPools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_lightingSets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };

        VkDescriptorSetLayout m_svgfTemporalLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_svgfTemporalPools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_svgfTemporalSets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkPipelineLayout m_svgfTemporalPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_svgfTemporalPipeline { VK_NULL_HANDLE };

        VkDescriptorSetLayout m_svgfATrousLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_svgfATrousPools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_svgfATrousSets[2][4]={};
        VkPipelineLayout m_svgfATrousPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_svgfATrousPipeline { VK_NULL_HANDLE };

        VkDescriptorSetLayout m_taaLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_taaPools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_taaSets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };

        VkDescriptorSetLayout m_blitLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_blitPools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_blitSets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };

        VkDescriptorSetLayout m_transparentSetLayout { VK_NULL_HANDLE };
        VkDescriptorPool m_transparentPool { VK_NULL_HANDLE };
        VkDescriptorSet m_transparentSet { VK_NULL_HANDLE };

        VkPipelineLayout m_graphicsPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_graphicsPipeline { VK_NULL_HANDLE };
        VkPipelineLayout m_transparentPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_transparentPipeline { VK_NULL_HANDLE };
        VkPipelineLayout m_lightingPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_lightingPipeline { VK_NULL_HANDLE };
        VkPipelineLayout m_taaPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_taaPipeline { VK_NULL_HANDLE };
        VkPipelineLayout m_blitPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_blitPipeline { VK_NULL_HANDLE };

        VkPipelineLayout m_debugLinePipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_debugLinePipeline { VK_NULL_HANDLE };
        std::unique_ptr<VulkanBuffer> m_debugLineBuffer;
        size_t m_maxDebugLineBufferSize = 2 * 1024 * 1024; // 2MB volume
        VkPipelineLayout m_editorGridPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_editorGridPipeline { VK_NULL_HANDLE };
        void CreateDebugLinePipeline();
        void CreateEditorGridPipeline();

        std::unordered_map<std::string, VulkanModelResource> m_modelCache;
        std::vector<std::unique_ptr<VulkanTexture>> m_globalTextures;
        std::vector<VkDescriptorImageInfo> m_globalImageInfos;
        std::unordered_map<std::string, uint32_t> m_globalTexturePathCache;
        VulkanFrameResource m_frameResource;

        struct OverrideMaterialBufferCacheEntry {
            std::string modelAssetPath;
            uint32_t componentRevision = 0;
            std::unique_ptr<VulkanBuffer> materialBuffer;
            uint64_t materialBufferAddress = 0;
        };
        std::unordered_map<uint32_t, OverrideMaterialBufferCacheEntry> m_overrideMaterialBuffers;

        std::unique_ptr<VulkanBuffer> m_rtInstanceBuffer;
        std::unique_ptr<VulkanBuffer> m_globalInstanceBuffer;
        std::unique_ptr<VulkanBuffer> m_pointLightBuffer;

        uint32_t m_width = 1280, m_height = 720;
        uint32_t m_frameIndex = 0;
        bool m_firstFrame = true, m_isFirstFrameRun = true;
        Matrix4x4 m_prevViewProj;
        std::unordered_map<uint32_t, Matrix4x4> m_prevModelMats;
        EditorViewportOverlayData m_editorOverlayData {};
        bool m_hasEditorOverlayData = false;
        RenderPipelinePreset m_renderPreset = RenderPipelinePreset::Stable;
        LightingProfile m_lightingProfile{};
        bool m_useManualLightingProfile = false;

        PFN_vkCmdDrawMeshTasksEXT m_CmdDrawMeshTasksEXT = nullptr;
        int m_viewX = 0, m_viewY = 0, m_viewW = 1280, m_viewH = 720;

        template <typename T>
        std::unique_ptr<VulkanBuffer> CreateBDABuffer(const std::vector<T>& data);

        void DestroyPipelines(VkDevice device);
        void DestroyDescriptors(VkDevice device);
        void DestroyAttachments();
        
        GBufferAttachment CreateAttachment(VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect);
        void CreateAttachments();
        void DestroyAttachment(GBufferAttachment& attachment);

        void BuildPipelines();
        void UpdateDescriptorSets();
        void UpdateLightingAccelerationStructureDescriptor(uint32_t ping, VkAccelerationStructureKHR tlasHandle);
        VulkanModelResource& GetOrLoadModel(const std::string& path);
        uint32_t GetOrLoadTextureIndex(const std::string& texturePath);
        uint64_t ResolveMaterialBufferAddress(Entity entity, const VulkanModelComponent& modelComp, const VulkanModelResource& modelResource);
        void ApplyMaterialOverride(MaterialData& targetMaterial, const VulkanMaterialSlotOverride& overrideData);
        void UpdateBindlessTextureDescriptor(uint32_t textureIndex);

        void TransitionImageLayout(VkCommandBuffer cmd, VkImage image, VkImageLayout oldLayout, VkImageLayout newLayout, VkImageAspectFlags aspectMask);
        float CreateHaltonSequence(uint32_t index, uint32_t base);
        LightingProfile ResolveLightingProfile(RenderPipelinePreset preset) const;
    };

} // namespace Lizeral
