#pragma once

#include "editor/overlay/EditorViewportOverlayTypes.h"
#include "runtime/function/ecs/registry.h"
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"
#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystemTypes.h"
#include "runtime/function/render/rhi/vulkan/VulkanBuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandPool.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandBuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanTLAS.h"
#include "runtime/function/render/passes/editor/EditorOverlayPasses.h"
#include "runtime/function/render/passes/mesh/MeshRenderPasses.h"
#include "runtime/function/render/passes/postprocess/PostProcessPasses.h"
#include "runtime/function/render/resources/RenderResourceCache.h"
#include "runtime/function/render/targets/RenderFrameTargets.h"
#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/core/math/matrix4.h"
#include "runtime/core/math/vector2.h"

#include <memory>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstdint>

namespace Lizeral {

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

        RenderFrameTargets m_frameTargets;

        DebugLinePass m_debugLinePass;
        EditorGridPass m_editorGridPass;
        EditorAxisPass m_editorAxisPass;
        TlasBuildPass m_tlasBuildPass;
        GBufferPass m_gBufferPass;
        TransparentPass m_transparentPass;
        LightingPass m_lightingPass;
        SvgfTemporalPass m_svgfTemporalPass;
        SvgfATrousPass m_svgfATrousPass;
        TaaPass m_taaPass;
        FinalBlitPass m_finalBlitPass;
        RenderResourceCache m_resourceCache;

        VulkanFrameResource m_frameResource;

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

        void BuildPipelines();
        void UpdateDescriptorSets();

        float CreateHaltonSequence(uint32_t index, uint32_t base);
        LightingProfile ResolveLightingProfile(RenderPipelinePreset preset) const;
    };

} // namespace Lizeral
