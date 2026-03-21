#pragma once

#include <vulkan/vulkan.h>
#include <vector>

namespace Lizeral {

    class VulkanDevice;

    class VulkanPipelineBuilder {
    public:
        VulkanPipelineBuilder();

        VulkanPipelineBuilder& AddShaderStage(VkShaderStageFlagBits stage, VkShaderModule shaderModule);
        VulkanPipelineBuilder& SetVertexInput(const VkPipelineVertexInputStateCreateInfo& vertexInputInfo);
        VulkanPipelineBuilder& SetInputAssembly(VkPrimitiveTopology topology, bool primitiveRestartEnable = false);
        VulkanPipelineBuilder& SetRasterization(VkPolygonMode polygonMode, VkCullModeFlags cullMode, VkFrontFace frontFace);
        VulkanPipelineBuilder& SetMultisampling(VkSampleCountFlagBits samples);
        VulkanPipelineBuilder& SetDepthStencil(bool depthTestEnable, bool depthWriteEnable, VkCompareOp depthCompareOp);
        
        VulkanPipelineBuilder& AddColorBlendAttachment(bool blendEnable, VkBlendOp colorBlendOp = VK_BLEND_OP_ADD, VkBlendFactor srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA, VkBlendFactor dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA);
        VulkanPipelineBuilder& DisableColorBlendAttachments(uint32_t count);

        VulkanPipelineBuilder& SetPipelineLayout(VkPipelineLayout pipelineLayout);

        VkPipeline Build(VulkanDevice* device, const std::vector<VkFormat>& colorFormats, VkFormat depthFormat);

    private:
        std::vector<VkPipelineShaderStageCreateInfo> m_shaderStages;
        VkPipelineVertexInputStateCreateInfo m_vertexInputInfo;
        VkPipelineInputAssemblyStateCreateInfo m_inputAssembly;
        VkPipelineViewportStateCreateInfo m_viewportState;
        VkPipelineRasterizationStateCreateInfo m_rasterizer;
        VkPipelineMultisampleStateCreateInfo m_multisampling;
        VkPipelineDepthStencilStateCreateInfo m_depthStencil;
        
        std::vector<VkPipelineColorBlendAttachmentState> m_colorBlendAttachments;
        VkPipelineColorBlendStateCreateInfo m_colorBlending;

        std::vector<VkDynamicState> m_dynamicStates;
        VkPipelineDynamicStateCreateInfo m_dynamicStateInfo;

        VkPipelineLayout m_pipelineLayout{ VK_NULL_HANDLE };
    };

} // namespace Lizeral