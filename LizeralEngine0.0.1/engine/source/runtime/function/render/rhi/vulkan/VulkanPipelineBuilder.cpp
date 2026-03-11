#include "VulkanPipelineBuilder.h"
#include "VulkanDevice.h"
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanPipelineBuilder::VulkanPipelineBuilder() {
        // 1. 只做最基础的 sType 初始化，清零内存
        m_vertexInputInfo = { VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO };
        m_inputAssembly = { VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO };
        m_viewportState = { VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO };
        m_rasterizer = { VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO };
        m_multisampling = { VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO };
        m_colorBlending = { VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO };
        m_depthStencil = { VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO };
        m_dynamicStateInfo = { VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO };

        // 2. 动态状态 (视口和裁剪永远是动态的，这是现代引擎的标配)
        m_viewportState.viewportCount = 1;
        m_viewportState.scissorCount = 1;
        m_dynamicStates = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };
        m_dynamicStateInfo.dynamicStateCount = static_cast<uint32_t>(m_dynamicStates.size());
        m_dynamicStateInfo.pDynamicStates = m_dynamicStates.data();
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::AddShaderStage(VkShaderStageFlagBits stage, VkShaderModule shaderModule) {
        VkPipelineShaderStageCreateInfo shaderStageInfo{ VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO };
        shaderStageInfo.stage = stage;
        shaderStageInfo.module = shaderModule;
        shaderStageInfo.pName = "main";
        m_shaderStages.push_back(shaderStageInfo);
        return *this;
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::SetVertexInput(const VkPipelineVertexInputStateCreateInfo& vertexInputInfo) {
        m_vertexInputInfo = vertexInputInfo;
        return *this;
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::SetInputAssembly(VkPrimitiveTopology topology, bool primitiveRestartEnable) {
        m_inputAssembly.topology = topology;
        m_inputAssembly.primitiveRestartEnable = primitiveRestartEnable ? VK_TRUE : VK_FALSE;
        return *this;
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::SetRasterization(VkPolygonMode polygonMode, VkCullModeFlags cullMode, VkFrontFace frontFace) {
        m_rasterizer.polygonMode = polygonMode;
        m_rasterizer.cullMode = cullMode;
        m_rasterizer.frontFace = frontFace;
        m_rasterizer.lineWidth = 1.0f;
        return *this;
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::SetMultisampling(VkSampleCountFlagBits samples) {
        m_multisampling.rasterizationSamples = samples;
        return *this;
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::SetDepthStencil(bool depthTestEnable, bool depthWriteEnable, VkCompareOp depthCompareOp) {
        m_depthStencil.depthTestEnable = depthTestEnable ? VK_TRUE : VK_FALSE;
        m_depthStencil.depthWriteEnable = depthWriteEnable ? VK_TRUE : VK_FALSE;
        m_depthStencil.depthCompareOp = depthCompareOp;
        return *this;
    }

    // ★ 逐个添加附件的混合模式
    VulkanPipelineBuilder& VulkanPipelineBuilder::AddColorBlendAttachment(bool blendEnable, VkBlendOp colorBlendOp, VkBlendFactor srcColorBlendFactor, VkBlendFactor dstColorBlendFactor) {
        VkPipelineColorBlendAttachmentState attachment{};
        attachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
        attachment.blendEnable = blendEnable ? VK_TRUE : VK_FALSE;
        if (blendEnable) {
            attachment.colorBlendOp = colorBlendOp;
            attachment.srcColorBlendFactor = srcColorBlendFactor;
            attachment.dstColorBlendFactor = dstColorBlendFactor;
            attachment.alphaBlendOp = VK_BLEND_OP_ADD;
            attachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
            attachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ZERO;
        }
        m_colorBlendAttachments.push_back(attachment);
        return *this;
    }

    // ★ 批量添加关闭混合的附件 (G-Buffer 极其需要，因为 G-Buffer 只是写数据，绝对不需要混合！)
    VulkanPipelineBuilder& VulkanPipelineBuilder::DisableColorBlendAttachments(uint32_t count) {
        for (uint32_t i = 0; i < count; i++) {
            AddColorBlendAttachment(false);
        }
        return *this;
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::SetPipelineLayout(VkPipelineLayout pipelineLayout) {
        m_pipelineLayout = pipelineLayout;
        return *this;
    }

    VkPipeline VulkanPipelineBuilder::Build(VulkanDevice* device, const std::vector<VkFormat>& colorFormats, VkFormat depthFormat) {
        if (m_pipelineLayout == VK_NULL_HANDLE) throw std::runtime_error("Cannot build pipeline: Pipeline Layout is null!");

        // 绑定颜色附件数组
        m_colorBlending.attachmentCount = static_cast<uint32_t>(m_colorBlendAttachments.size());
        m_colorBlending.pAttachments = m_colorBlendAttachments.data();

        // ★ 核心改变：接收多个 Color Format
        VkPipelineRenderingCreateInfo renderingInfo{ VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO };
        renderingInfo.colorAttachmentCount = static_cast<uint32_t>(colorFormats.size());
        renderingInfo.pColorAttachmentFormats = colorFormats.data();
        renderingInfo.depthAttachmentFormat = depthFormat;

        VkGraphicsPipelineCreateInfo pipelineInfo{ VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO };
        pipelineInfo.pNext = &renderingInfo; 
        pipelineInfo.stageCount = static_cast<uint32_t>(m_shaderStages.size());
        pipelineInfo.pStages = m_shaderStages.data(); 

        pipelineInfo.pVertexInputState = &m_vertexInputInfo;
        pipelineInfo.pInputAssemblyState = &m_inputAssembly;
        pipelineInfo.pViewportState = &m_viewportState;
        pipelineInfo.pRasterizationState = &m_rasterizer;
        pipelineInfo.pMultisampleState = &m_multisampling;
        pipelineInfo.pDepthStencilState = &m_depthStencil;
        pipelineInfo.pColorBlendState = &m_colorBlending;
        pipelineInfo.pDynamicState = &m_dynamicStateInfo; 

        pipelineInfo.layout = m_pipelineLayout;
        pipelineInfo.renderPass = VK_NULL_HANDLE; 

        VkPipeline graphicsPipeline = VK_NULL_HANDLE;
        if (vkCreateGraphicsPipelines(device->GetNativeDevice(), VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &graphicsPipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create graphics pipeline!");
        }

        std::cout << "[VulkanPipeline] Graphics Pipeline built successfully for Dynamic Rendering!" << std::endl;
        return graphicsPipeline;
    }

} // namespace Lizeral