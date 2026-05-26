#pragma once

#include "editor/overlay/EditorViewportOverlayTypes.h"
#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/function/render/RenderFrameContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanBuffer.h"

#include <memory>
#include <vector>
#include <vulkan/vulkan.h>

namespace Lizeral {

    class VulkanDevice;

    class DebugLinePass {
    public:
        void Initialize(VulkanDevice* device, VkFormat colorFormat, VkDescriptorSetLayout sceneDepthSetLayout);
        void Shutdown(VkDevice device);
        void Render(const RenderFrameContext& ctx, const std::vector<DebugLineVertex>& lines);

    private:
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
        std::unique_ptr<VulkanBuffer> m_lineBuffer;
        size_t m_maxLineBufferSize = 2 * 1024 * 1024;
    };

    class EditorGridPass {
    public:
        void Initialize(VulkanDevice* device, VkFormat colorFormat, VkDescriptorSetLayout sceneDepthSetLayout);
        void Shutdown(VkDevice device);
        void Render(const RenderFrameContext& ctx, const EditorViewportOverlayData& overlay);

    private:
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

    class EditorAxisPass {
    public:
        void Initialize(VulkanDevice* device, VkFormat colorFormat, VkDescriptorSetLayout sceneDepthSetLayout);
        void Shutdown(VkDevice device);
        void Render(const RenderFrameContext& ctx, const EditorViewportOverlayData& overlay);

    private:
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
    };

} // namespace Lizeral
