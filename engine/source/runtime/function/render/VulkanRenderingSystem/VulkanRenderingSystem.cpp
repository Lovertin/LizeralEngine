#include "VulkanRenderingSystem.h"
#include "runtime/function/render/draw/RenderDrawPacket.h"
#include "runtime/function/render/passes/RenderPassUtils.h"
#include "runtime/function/render/scene/RenderSceneSnapshot.h"

#include <iostream>
#include <cmath>

namespace Lizeral {

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

    float VulkanRenderingSystem::CreateHaltonSequence(uint32_t index, uint32_t base) {
        float f = 1.0f; float result = 0.0f; uint32_t i = index;
        while (i > 0) {
            f = f / static_cast<float>(base);
            result = result + f * static_cast<float>(i % base);
            i = static_cast<uint32_t>(std::floor(static_cast<float>(i) / static_cast<float>(base)));
        }
        return result;
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
            m_resourceCache.Initialize(m_device, m_resourceCommandPool.get());
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

        m_frameTargets.Initialize(m_device, m_width, m_height);

        BuildPipelines();

        m_firstFrame = true;
        m_isFirstFrameRun = true;

        std::cout << "[VulkanRenderingSystem] Initialization Complete." << std::endl;
    }

    void VulkanRenderingSystem::Shutdown() {
        if (!m_device) return;
        vkDeviceWaitIdle(m_device->GetNativeDevice());
        VkDevice device = m_device->GetNativeDevice();

        std::cout << "[RenderingSystem] Shutting down..." << std::endl;

        m_resourceCache.Clear();

        if (m_rtInstanceBuffer) m_rtInstanceBuffer.reset();
        if (m_globalInstanceBuffer) m_globalInstanceBuffer.reset(); 
        if (m_frameResource.buffer) m_frameResource.buffer.reset(); 
        if (m_pointLightBuffer) m_pointLightBuffer.reset();
        m_debugLinePass.Shutdown(device);
        m_editorGridPass.Shutdown(device);
        m_editorAxisPass.Shutdown(device);

        DestroyPipelines(device);
        DestroyDescriptors(device);
        m_frameTargets.Shutdown();
        if (m_tlas) m_tlas.reset();
        if (m_resourceCommandPool) m_resourceCommandPool.reset();

        std::cout << "[RenderingSystem] Shutdown complete." << std::endl;
    }

    void VulkanRenderingSystem::DestroyPipelines(VkDevice device){
        m_gBufferPass.Shutdown(device);
        m_transparentPass.Shutdown(device);
        m_lightingPass.Shutdown(device);
        m_svgfTemporalPass.Shutdown(device);
        m_svgfATrousPass.Shutdown(device);
        m_taaPass.Shutdown(device);
        m_finalBlitPass.Shutdown(device);
        m_debugLinePass.Shutdown(device);
        m_editorGridPass.Shutdown(device);
        m_editorAxisPass.Shutdown(device);
    }

    void VulkanRenderingSystem::DestroyDescriptors(VkDevice device){
        (void)device;
    }

    void VulkanRenderingSystem::Resize(uint32_t width, uint32_t height) {
        if (width == 0 || height == 0) return;
        if (width == m_width && height == m_height) return;
        vkDeviceWaitIdle(m_device->GetNativeDevice());

        m_width = width;
        m_height = height;

        m_frameTargets.Resize(width, height);

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

    void VulkanRenderingSystem::BuildPipelines() {
        std::cout << "[RenderingSystem] Building Pipelines..." << std::endl;
        m_gBufferPass.Initialize(m_device, m_resourceCache, m_frameTargets);
        m_transparentPass.Initialize(m_device, m_renderer->GetSwapchainFormat(), m_gBufferPass.GetBindlessSetLayout());
        LightingPassSettings lightingSettings{};
        lightingSettings.giQualityLevel = m_lightingProfile.giQualityLevel;
        lightingSettings.shadowQualityLevel = m_lightingProfile.shadowQualityLevel;
        m_lightingPass.Initialize(m_device, m_gBufferPass.GetBindlessSetLayout(), lightingSettings);
        m_svgfTemporalPass.Initialize(m_device);
        m_svgfATrousPass.Initialize(m_device);
        m_taaPass.Initialize(m_device);
        m_finalBlitPass.Initialize(m_device, m_renderer->GetSwapchainFormat());

        UpdateDescriptorSets();

        VkFormat swapchainFormat = m_renderer->GetSwapchainFormat();
        m_editorGridPass.Initialize(m_device, swapchainFormat, m_transparentPass.GetSceneDepthSetLayout());
        m_editorAxisPass.Initialize(m_device, swapchainFormat, m_transparentPass.GetSceneDepthSetLayout());
        m_debugLinePass.Initialize(m_device, swapchainFormat, m_transparentPass.GetSceneDepthSetLayout());
        std::cout << "[RenderingSystem] Pipelines Built Successfully." << std::endl;
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

        RenderSceneSnapshot scene = ExtractRenderSceneSnapshot(registry, m_width, m_height, ndcJitterX, ndcJitterY);

        Matrix4x4 viewMat = scene.camera.view;
        Matrix4x4 projMatUnjittered = scene.camera.projection;
        Matrix4x4 projMatJittered = scene.camera.jitteredProjection;
        Vector3 cameraPos = scene.camera.position;
        Vector3 cameraRight = scene.camera.right;
        Vector3 cameraUp = scene.camera.up;
        Vector3 cameraForward = scene.camera.forward;
        const Vector3& lightDir = scene.directionalLight.direction;
        const Vector3& lightColor = scene.directionalLight.color;
        const float lightIntensity = scene.directionalLight.intensity;

        Matrix4x4 currentVP = projMatUnjittered * viewMat; //projMatUnjittered is important !!!
        if (m_isFirstFrameRun) {
            m_prevViewProj = currentVP;
        }

        RenderFrameContext frameCtx{};
        frameCtx.cmd = cmd;
        frameCtx.viewport = viewport;
        frameCtx.scissor = scissor;
        frameCtx.width = m_width;
        frameCtx.height = m_height;
        frameCtx.frameIndex = m_frameIndex;
        frameCtx.view = viewMat;
        frameCtx.projection = projMatUnjittered;
        frameCtx.viewProj = currentVP;
        frameCtx.cameraPos = cameraPos;
        frameCtx.cameraRight = cameraRight;
        frameCtx.cameraUp = cameraUp;
        frameCtx.cameraForward = cameraForward;
        frameCtx.sceneDepthSet = m_transparentPass.GetSceneDepthSet();

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

        RenderDrawPacket drawPacket = BuildRenderDrawPacket(scene, m_resourceCache);
        if (!drawPacket.rtInstanceDescs.empty()) {
            m_rtInstanceBuffer->WriteData(
                drawPacket.rtInstanceDescs.data(),
                drawPacket.rtInstanceDescs.size() * sizeof(RTInstanceDesc)
            );
        }

        m_tlasBuildPass.Render(frameCtx, *m_tlas, ping, drawPacket);

        // hot swap Descriptor
        VkAccelerationStructureKHR currentTlas = m_tlas->GetHandle(ping);
        m_lightingPass.UpdateAccelerationStructureDescriptor(m_device, ping, currentTlas);

        //  G-Buffer Pass
        uint64_t baseInstanceAddr = m_globalInstanceBuffer->GetDeviceAddress();
        WriteRenderDrawPacketInstanceData(drawPacket, m_isFirstFrameRun, m_prevModelMats, *m_globalInstanceBuffer);
        m_gBufferPass.Render(frameCtx, m_frameTargets, drawPacket, m_frameResource.addr, baseInstanceAddr, m_firstFrame, m_CmdDrawMeshTasksEXT);

        uint32_t lightCount = static_cast<uint32_t>(scene.pointLights.size());
        uint32_t maxLightCapacity = 1024; 
        if (lightCount > maxLightCapacity) {
            std::cerr << "[Warning] Point lights count exceeds capacity (1024), truncating!" << std::endl;
            lightCount = maxLightCapacity;
        }
        if (lightCount > 0) { 
            m_pointLightBuffer->WriteData(scene.pointLights.data(), lightCount * sizeof(RenderPointLightSnapshot));
        }

        m_prevViewProj = currentVP;
        if (m_isFirstFrameRun) {
            TransitionImageLayout(cmd, m_frameTargets.directLight.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.noisyGI.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.denoisedGI.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.denoisedGITemp.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.history[0].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.history[1].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.giHistory[0].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.giHistory[1].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.momentsHistory[0].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            TransitionImageLayout(cmd, m_frameTargets.momentsHistory[1].image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        }

        FullscreenPassPushConstants lightPc{};
        lightPc.frameDataAddr = m_frameResource.addr;
        lightPc.instanceDescAddr = m_rtInstanceBuffer->GetDeviceAddress();
        lightPc.pointLightsAddr = m_pointLightBuffer ? m_pointLightBuffer->GetDeviceAddress() : 0;
        lightPc.pointLightCount = lightCount;


        // 1. Lighting Pass
        m_lightingPass.Render(frameCtx, m_frameTargets, ping, m_gBufferPass.GetBindlessSet(), lightPc);

        // 2. SVGF Temporal Pass
        m_svgfTemporalPass.Render(frameCtx, m_frameTargets, ping, lightPc);

        // 3. SVGF A-Trous Pass
        if (m_isFirstFrameRun) {
            TransitionImageLayout(cmd, m_frameTargets.denoisedGITemp.image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        }

        m_svgfATrousPass.Render(frameCtx, m_frameTargets, ping, lightPc);

        // 4. TAA Pass
        m_taaPass.Render(frameCtx, m_frameTargets, ping, lightPc);

        // 5. Blit Pass
        m_renderer->BeginRendering(cmd); 
        vkCmdSetViewport(cmd, 0, 1, &viewport); vkCmdSetScissor(cmd, 0, 1, &scissor);
        m_finalBlitPass.Render(frameCtx, ping, lightPc);

        if (activeOverlay) {
            m_editorGridPass.Render(frameCtx, *activeOverlay);
        }

        m_transparentPass.Render(
            frameCtx,
            drawPacket,
            m_gBufferPass.GetBindlessSet(),
            m_frameResource.addr,
            baseInstanceAddr,
            m_CmdDrawMeshTasksEXT
        );

        // 6. Editor overlay lines and viewport gizmo.
        m_debugLinePass.Render(frameCtx, activeDebugLines);
        if (activeOverlay) {
            m_editorAxisPass.Render(frameCtx, *activeOverlay);
        }

        m_renderer->EndRendering(cmd);

        m_renderer->EndFrame();
        m_isFirstFrameRun = m_firstFrame = false;
    }


    void VulkanRenderingSystem::UpdateDescriptorSets() {
        m_transparentPass.UpdateDescriptors(m_device, m_frameTargets);
        for (uint32_t frameSlot = 0; frameSlot < 2; ++frameSlot) {
            VkAccelerationStructureKHR tlasHandle = m_tlas ? m_tlas->GetHandle(frameSlot) : VK_NULL_HANDLE;
            m_lightingPass.UpdateDescriptors(m_device, m_frameTargets, frameSlot, tlasHandle);
            m_svgfTemporalPass.UpdateDescriptors(m_device, m_frameTargets, frameSlot);
            m_svgfATrousPass.UpdateDescriptors(m_device, m_frameTargets, frameSlot);
            m_taaPass.UpdateDescriptors(m_device, m_frameTargets, frameSlot);
            m_finalBlitPass.UpdateDescriptors(m_device, m_frameTargets, frameSlot);
        }
    }

} // namespace Lizeral
