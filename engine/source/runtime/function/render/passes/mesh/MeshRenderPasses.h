#pragma once

#include "runtime/function/render/RenderFrameContext.h"

#include <vulkan/vulkan.h>

namespace Lizeral {

    class RenderFrameTargets;
    class RenderResourceCache;
    class VulkanDevice;
    class VulkanTLAS;
    struct RenderDrawPacket;

    class TlasBuildPass {
    public:
        void Render(const RenderFrameContext& ctx, VulkanTLAS& tlas, uint32_t frameSlot, const RenderDrawPacket& packet);
    };

    class GBufferPass {
    public:
        void Initialize(VulkanDevice* device, RenderResourceCache& resourceCache, const RenderFrameTargets& targets);
        void Shutdown(VkDevice device);

        VkDescriptorSetLayout GetBindlessSetLayout() const { return m_bindlessSetLayout; }
        VkDescriptorSet GetBindlessSet() const { return m_bindlessSet; }

        void Render(
            const RenderFrameContext& ctx,
            RenderFrameTargets& targets,
            const RenderDrawPacket& packet,
            uint64_t frameDataAddr,
            uint64_t baseInstanceAddr,
            bool firstFrame,
            PFN_vkCmdDrawMeshTasksEXT cmdDrawMeshTasks
        );

    private:
        VkDescriptorSetLayout m_bindlessSetLayout { VK_NULL_HANDLE };
        VkDescriptorPool m_bindlessPool { VK_NULL_HANDLE };
        VkDescriptorSet m_bindlessSet { VK_NULL_HANDLE };
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

    class TransparentPass {
    public:
        void Initialize(VulkanDevice* device, VkFormat colorFormat, VkDescriptorSetLayout bindlessSetLayout);
        void Shutdown(VkDevice device);
        void UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets);

        VkDescriptorSetLayout GetSceneDepthSetLayout() const { return m_sceneDepthSetLayout; }
        VkDescriptorSet GetSceneDepthSet() const { return m_sceneDepthSet; }

        void Render(
            const RenderFrameContext& ctx,
            const RenderDrawPacket& packet,
            VkDescriptorSet bindlessSet,
            uint64_t frameDataAddr,
            uint64_t baseInstanceAddr,
            PFN_vkCmdDrawMeshTasksEXT cmdDrawMeshTasks
        );

    private:
        VkDescriptorSetLayout m_sceneDepthSetLayout { VK_NULL_HANDLE };
        VkDescriptorPool m_sceneDepthPool { VK_NULL_HANDLE };
        VkDescriptorSet m_sceneDepthSet { VK_NULL_HANDLE };
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

} // namespace Lizeral
