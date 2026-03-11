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
        
        // ★ 核心改变：支持添加多个颜色混合附件 (为 MRT 准备)
        VulkanPipelineBuilder& AddColorBlendAttachment(bool blendEnable, VkBlendOp colorBlendOp = VK_BLEND_OP_ADD, VkBlendFactor srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA, VkBlendFactor dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA);
        VulkanPipelineBuilder& DisableColorBlendAttachments(uint32_t count); // 快速添加多个关闭混合的附件 (G-Buffer 最常用)

        VulkanPipelineBuilder& SetPipelineLayout(VkPipelineLayout pipelineLayout);

        // ★ 修改 Build 函数签名，接收一个颜色格式的数组！
        VkPipeline Build(VulkanDevice* device, const std::vector<VkFormat>& colorFormats, VkFormat depthFormat);

    private:
        std::vector<VkPipelineShaderStageCreateInfo> m_shaderStages;
        VkPipelineVertexInputStateCreateInfo m_vertexInputInfo;
        VkPipelineInputAssemblyStateCreateInfo m_inputAssembly;
        VkPipelineViewportStateCreateInfo m_viewportState;
        VkPipelineRasterizationStateCreateInfo m_rasterizer;
        VkPipelineMultisampleStateCreateInfo m_multisampling;
        VkPipelineDepthStencilStateCreateInfo m_depthStencil;
        
        // ★ 核心改变：变成数组
        std::vector<VkPipelineColorBlendAttachmentState> m_colorBlendAttachments;
        VkPipelineColorBlendStateCreateInfo m_colorBlending;

        std::vector<VkDynamicState> m_dynamicStates;
        VkPipelineDynamicStateCreateInfo m_dynamicStateInfo;

        VkPipelineLayout m_pipelineLayout{ VK_NULL_HANDLE };
    };

} // namespace Lizeral