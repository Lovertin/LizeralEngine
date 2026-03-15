#pragma once

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
#include "runtime/core/math/matrix4.h"

#include <memory>
#include <unordered_map>
#include <vector>
#include <string>

namespace Lizeral {

    // 从 sandbox 迁移过来的核心数据结构
    struct GBufferAttachment {
        VkImage image = VK_NULL_HANDLE;
        VmaAllocation allocation = VK_NULL_HANDLE;
        VkImageView view = VK_NULL_HANDLE;
        VkFormat format;
    };

    class VulkanRenderingSystem {
    public:
        VulkanRenderingSystem() = default;
        ~VulkanRenderingSystem() {  }

        // ★ 初始化：由引擎核心在创建好 Device 和 Renderer 后注入
        void Initialize(VulkanContext* context, VulkanDevice* device, VulkanRenderer* renderer, uint32_t width, uint32_t height);

        // ★ 每帧核心渲染循环
        void Tick(Registry& registry, float deltaTime);

        // ★ 清理资源
        void Shutdown();

        // 供外部编辑器或窗口系统在拉伸窗口时调用
        void Resize(uint32_t width, uint32_t height);

        void SetViewport(int x, int y, int w, int h) {
            m_viewX = x; m_viewY = y; m_viewW = w; m_viewH = h;
        }
        // Vulkan 中通常不需要手动 SetDefaultFBO，留空兼容旧代码接口即可
        void SetDefaultFBO(unsigned int fbo) { /* VulkanRenderer handles target */ }

    private:
        // --- 核心 RHI 句柄 (不拥有生命周期，由外部传入) ---
        VulkanContext* m_context { nullptr };
        VulkanDevice* m_device { nullptr };
        VulkanRenderer* m_renderer { nullptr };

        // --- 全局基础设施 ---
        std::unique_ptr<VulkanCommandPool> m_resourceCommandPool;
        std::unique_ptr<VulkanTLAS> m_tlas;

        // --- G-Buffer 附件 ---
        GBufferAttachment m_gAlbedoMetallic;
        GBufferAttachment m_gNormalRoughness;
        GBufferAttachment m_gDepth;
        GBufferAttachment m_gVelocity;
        GBufferAttachment m_gDirectLight;
        GBufferAttachment m_gNoisyGI;
        GBufferAttachment m_gDenoisedGI;
        GBufferAttachment m_gHistory[2];
        VkSampler m_gBufferSampler { VK_NULL_HANDLE };

        // --- 描述符与管线资源 ---
        // Bindless (Mesh)
        VkDescriptorSetLayout m_descriptorSetLayout { VK_NULL_HANDLE };
        VkDescriptorPool m_descriptorPool { VK_NULL_HANDLE };
        VkDescriptorSet m_descriptorSet { VK_NULL_HANDLE };

        // Lighting Pass (Ping-pong)
        VkDescriptorSetLayout m_lightingLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_lightPools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_lightingSets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };

        // Denoise Pass (Ping-pong)
        VkDescriptorSetLayout m_denoiseLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_denoisePools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_denoiseSets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };

        // TAA Pass (Ping-pong)
        VkDescriptorSetLayout m_taaLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_taaPools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_taaSets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };

        // Blit Pass (Ping-pong)
        VkDescriptorSetLayout m_blitLayouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_blitPools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_blitSets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };

        // 管线布局与管线对象
        VkPipelineLayout m_graphicsPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_graphicsPipeline { VK_NULL_HANDLE };

        VkPipelineLayout m_lightingPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_lightingPipeline { VK_NULL_HANDLE };

        VkPipelineLayout m_denoisePipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_denoisePipeline { VK_NULL_HANDLE };

        VkPipelineLayout m_taaPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_taaPipeline { VK_NULL_HANDLE };

        VkPipelineLayout m_blitPipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_blitPipeline { VK_NULL_HANDLE };

        // --- 资源管理台账 ---
        std::unordered_map<std::string, VulkanModelResource> m_modelCache;
        std::vector<std::unique_ptr<VulkanTexture>> m_globalTextures;
        std::vector<VkDescriptorImageInfo> m_globalImageInfos;

        // 动态光追台账 Buffer
        std::unique_ptr<VulkanBuffer> m_rtInstanceBuffer;

        // --- 状态数据 ---
        uint32_t m_width = 1280;
        uint32_t m_height = 720;
        uint32_t m_frameIndex = 0;
        bool m_firstFrame = true;

        // --- 渲染状态缓存 (TAA 与 运动向量需要) ---
        Matrix4x4 m_prevViewProj;
        std::unordered_map<uint32_t, Matrix4x4> m_prevModelMats; // Entity ID -> Matrix
        bool m_isFirstFrameRun = true;

        // --- 动态扩展函数指针 ---
        PFN_vkCmdDrawMeshTasksEXT m_CmdDrawMeshTasksEXT = nullptr;

        // --- 视口与裁剪区域 (由 Editor 传入) ---
        int m_viewX = 0, m_viewY = 0, m_viewW = 1280, m_viewH = 720;

        // --- 内部辅助方法 ---
        template <typename T>
        std::unique_ptr<VulkanBuffer> CreateBDABuffer(const std::vector<T>& data);
        
        GBufferAttachment CreateAttachment(VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect);
        void DestroyAttachment(GBufferAttachment& attachment);

        void BuildPipelines();
        VulkanModelResource& GetOrLoadModel(const std::string& path);

        void TransitionImageLayout(VkCommandBuffer cmd, VkImage image, VkImageLayout oldLayout, VkImageLayout newLayout, VkImageAspectFlags aspectMask);
        float CreateHaltonSequence(uint32_t index, uint32_t base);
    };

} // namespace Lizeral