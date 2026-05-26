#pragma once

#include "runtime/function/render/RenderFrameContext.h"
#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystemTypes.h"
#include "runtime/function/render/targets/RenderFrameTargets.h"

#include <vulkan/vulkan.h>

namespace Lizeral {

    class VulkanDevice;

    struct LightingPassSettings {
        int32_t giQualityLevel = 0;
        int32_t shadowQualityLevel = 0;
    };

    class LightingPass {
    public:
        void Initialize(VulkanDevice* device, VkDescriptorSetLayout bindlessSetLayout, const LightingPassSettings& settings);
        void Shutdown(VkDevice device);
        void UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot, VkAccelerationStructureKHR tlasHandle);
        void UpdateAccelerationStructureDescriptor(VulkanDevice* device, uint32_t frameSlot, VkAccelerationStructureKHR tlasHandle);
        void Render(const RenderFrameContext& ctx, RenderFrameTargets& targets, uint32_t frameSlot, VkDescriptorSet bindlessSet, const FullscreenPassPushConstants& pushConstants);

    private:
        VkDescriptorSetLayout m_layouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_pools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_sets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

    class SvgfTemporalPass {
    public:
        void Initialize(VulkanDevice* device);
        void Shutdown(VkDevice device);
        void UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot);
        void Render(const RenderFrameContext& ctx, RenderFrameTargets& targets, uint32_t frameSlot, const FullscreenPassPushConstants& pushConstants);

    private:
        VkDescriptorSetLayout m_layouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_pools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_sets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

    class SvgfATrousPass {
    public:
        void Initialize(VulkanDevice* device);
        void Shutdown(VkDevice device);
        void UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot);
        void Render(const RenderFrameContext& ctx, RenderFrameTargets& targets, uint32_t frameSlot, FullscreenPassPushConstants pushConstants);

    private:
        VkDescriptorSetLayout m_layouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_pools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_sets[2][4] {};
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

    class TaaPass {
    public:
        void Initialize(VulkanDevice* device);
        void Shutdown(VkDevice device);
        void UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot);
        void Render(const RenderFrameContext& ctx, RenderFrameTargets& targets, uint32_t frameSlot, const FullscreenPassPushConstants& pushConstants);

    private:
        VkDescriptorSetLayout m_layouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_pools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_sets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

    class FinalBlitPass {
    public:
        void Initialize(VulkanDevice* device, VkFormat colorFormat);
        void Shutdown(VkDevice device);
        void UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot);
        void Render(const RenderFrameContext& ctx, uint32_t frameSlot, const FullscreenPassPushConstants& pushConstants);

    private:
        VkDescriptorSetLayout m_layouts[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorPool m_pools[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkDescriptorSet m_sets[2] { VK_NULL_HANDLE, VK_NULL_HANDLE };
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

} // namespace Lizeral
