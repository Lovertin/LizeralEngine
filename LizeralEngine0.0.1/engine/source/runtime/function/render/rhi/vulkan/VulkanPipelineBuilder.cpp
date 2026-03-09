#include "VulkanPipelineBuilder.h"
#include "VulkanDevice.h"
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanPipelineBuilder::VulkanPipelineBuilder() {
        // =========================================================
        // 赋予现代引擎最常用的默认值 (所有没被用户显式 Set 的状态)
        // =========================================================

        // 1. 默认没有任何顶点输入 (完美契合 BDA 和 Mesh Shader)
        m_vertexInputInfo = {};
        m_vertexInputInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;

        // 2. 默认画三角形
        m_inputAssembly = {};
        m_inputAssembly.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
        m_inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
        m_inputAssembly.primitiveRestartEnable = VK_FALSE;

        // 3. ★ 核心现代特性：永远开启动态视口和裁剪
        m_viewportState = {};
        m_viewportState.sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO;
        m_viewportState.viewportCount = 1; // 数量必须填，但指针留空给动态状态
        m_viewportState.scissorCount = 1;

        m_dynamicStates = {
            VK_DYNAMIC_STATE_VIEWPORT,
            VK_DYNAMIC_STATE_SCISSOR
        };
        m_dynamicStateInfo = {};
        m_dynamicStateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO;
        m_dynamicStateInfo.dynamicStateCount = static_cast<uint32_t>(m_dynamicStates.size());
        m_dynamicStateInfo.pDynamicStates = m_dynamicStates.data();

        // 4. 默认光栅化：填充模式，不剔除背面 (方便我们第一次测三角形)
        m_rasterizer = {};
        m_rasterizer.sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO;
        m_rasterizer.depthClampEnable = VK_FALSE;
        m_rasterizer.rasterizerDiscardEnable = VK_FALSE;
        m_rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
        m_rasterizer.lineWidth = 1.0f;
        m_rasterizer.cullMode = VK_CULL_MODE_NONE; // 暂时关闭剔除
        m_rasterizer.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE; // 逆时针为正面

        // 5. 默认抗锯齿：关闭 (1 个采样点)
        m_multisampling = {};
        m_multisampling.sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO;
        m_multisampling.sampleShadingEnable = VK_FALSE;
        m_multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

        // 6. 默认颜色混合：覆盖 (不透明)
        m_colorBlendAttachment = {};
        m_colorBlendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
        m_colorBlendAttachment.blendEnable = VK_FALSE; // 不开启混合

        m_colorBlending = {};
        m_colorBlending.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
        m_colorBlending.logicOpEnable = VK_FALSE;
        m_colorBlending.attachmentCount = 1;
        m_colorBlending.pAttachments = &m_colorBlendAttachment;

        // 7. 默认深度测试：关闭 (画基本 UI 或三角形时不需要)
        m_depthStencil = {};
        m_depthStencil.sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO;
        m_depthStencil.depthTestEnable = VK_TRUE;           // 开启深度测试
        m_depthStencil.depthWriteEnable = VK_TRUE;          // 允许写入深度缓冲（极其重要，否则前面的挡不住后面的）
        m_depthStencil.depthCompareOp = VK_COMPARE_OP_LESS; // 深度值越小（越近）越优先显示
        m_depthStencil.depthBoundsTestEnable = VK_FALSE;
        m_depthStencil.stencilTestEnable = VK_FALSE;
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::AddShaderStage(VkShaderStageFlagBits stage, VkShaderModule shaderModule) {
        VkPipelineShaderStageCreateInfo shaderStageInfo{};
        shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        shaderStageInfo.stage = stage;
        shaderStageInfo.module = shaderModule;
        shaderStageInfo.pName = "main"; // 对应 glsl 里的入口函数名

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

    VulkanPipelineBuilder& VulkanPipelineBuilder::SetColorBlend(bool blendEnable, VkBlendOp colorBlendOp, VkBlendFactor srcColorBlendFactor, VkBlendFactor dstColorBlendFactor) {
        m_colorBlendAttachment.blendEnable = blendEnable ? VK_TRUE : VK_FALSE;
        if (blendEnable) {
            m_colorBlendAttachment.colorBlendOp = colorBlendOp;
            m_colorBlendAttachment.srcColorBlendFactor = srcColorBlendFactor;
            m_colorBlendAttachment.dstColorBlendFactor = dstColorBlendFactor;
            // Alpha通道通常这样混合即可
            m_colorBlendAttachment.alphaBlendOp = VK_BLEND_OP_ADD;
            m_colorBlendAttachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
            m_colorBlendAttachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ZERO;
        }
        return *this;
    }

    VulkanPipelineBuilder& VulkanPipelineBuilder::SetPipelineLayout(VkPipelineLayout pipelineLayout) {
        m_pipelineLayout = pipelineLayout;
        return *this;
    }

    VkPipeline VulkanPipelineBuilder::Build(VulkanDevice* device, VkFormat colorFormat, VkFormat depthFormat) {
        if (m_pipelineLayout == VK_NULL_HANDLE) {
            throw std::runtime_error("Cannot build pipeline: Pipeline Layout is null!");
        }

        // =========================================================
        // ★ 动态渲染核心魔法：告诉管线我们要画什么格式的图片
        // =========================================================
        VkPipelineRenderingCreateInfo renderingInfo{};
        renderingInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO;
        renderingInfo.colorAttachmentCount = 1;
        renderingInfo.pColorAttachmentFormats = &colorFormat;
        renderingInfo.depthAttachmentFormat = depthFormat;

        // 最终组装装配线！
        VkGraphicsPipelineCreateInfo pipelineInfo{};
        pipelineInfo.sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO;
        pipelineInfo.pNext = &renderingInfo; // ★ 挂载动态渲染格式信息！

        // ★ 就是这里！刚刚因为缩略写法不小心丢失的代码：
        pipelineInfo.stageCount = static_cast<uint32_t>(m_shaderStages.size());
        pipelineInfo.pStages = m_shaderStages.data();  // 把 Shader 数组传给管线！

        // 其他状态机装配
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
        pipelineInfo.subpass = 0;             

        VkPipeline graphicsPipeline = VK_NULL_HANDLE;
        if (vkCreateGraphicsPipelines(device->GetNativeDevice(), VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &graphicsPipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create graphics pipeline!");
        }

        std::cout << "[VulkanPipeline] Graphics Pipeline built successfully for Dynamic Rendering!" << std::endl;
        return graphicsPipeline;
    }

} // namespace Lizeral