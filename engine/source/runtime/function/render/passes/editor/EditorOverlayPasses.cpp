#include "runtime/function/render/passes/editor/EditorOverlayPasses.h"

#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/rhi/vulkan/VulkanShaderUtils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace Lizeral {
namespace {

    struct EditorGridPushConstants {
        Matrix4x4 viewProj;
        Matrix4x4 invViewProj;
        Vector4 cameraPosAndPlaneHeight;
        Vector4 viewportSizeAndSpacing;
        Vector4 fadeAndOpacity;
    };

    struct EditorAxisPushConstants {
        Matrix4x4 viewProj;
        Vector4 viewportSizeAndThickness;
        Vector4 axisStart;
        Vector4 axisEnd;
        Vector4 axisColor;
    };

    float SmoothStep(float edge0, float edge1, float value) {
        float t = std::clamp((value - edge0) / (edge1 - edge0), 0.0f, 1.0f);
        return t * t * (3.0f - 2.0f * t);
    }

} // namespace

    void DebugLinePass::Initialize(VulkanDevice* device, VkFormat colorFormat, VkDescriptorSetLayout sceneDepthSetLayout) {
        VkDevice nativeDevice = device->GetNativeDevice();
        const std::string shaderDir = LIZERAL_SHADER_DIR;

        VkShaderModule vertShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "debug_line_vert.spv"));
        VkShaderModule fragShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "debug_line_frag.spv"));

        VkPipelineShaderStageCreateInfo shaderStages[] = {
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_VERTEX_BIT, vertShader, "main", nullptr},
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_FRAGMENT_BIT, fragShader, "main", nullptr}
        };

        VkVertexInputBindingDescription bindingDesc{};
        bindingDesc.binding = 0;
        bindingDesc.stride = sizeof(DebugLineVertex);
        bindingDesc.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;

        VkVertexInputAttributeDescription attrDesc[2] = {};
        attrDesc[0].binding = 0;
        attrDesc[0].location = 0;
        attrDesc[0].format = VK_FORMAT_R32G32B32_SFLOAT;
        attrDesc[0].offset = offsetof(DebugLineVertex, position);
        attrDesc[1].binding = 0;
        attrDesc[1].location = 1;
        attrDesc[1].format = VK_FORMAT_R32G32B32_SFLOAT;
        attrDesc[1].offset = offsetof(DebugLineVertex, color);

        VkPipelineVertexInputStateCreateInfo vertexInputInfo{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
        vertexInputInfo.vertexBindingDescriptionCount = 1;
        vertexInputInfo.pVertexBindingDescriptions = &bindingDesc;
        vertexInputInfo.vertexAttributeDescriptionCount = 2;
        vertexInputInfo.pVertexAttributeDescriptions = attrDesc;

        VkPipelineInputAssemblyStateCreateInfo inputAssembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
        inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_LINE_LIST;

        VkPipelineViewportStateCreateInfo viewportState{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
        viewportState.viewportCount = 1;
        viewportState.scissorCount = 1;

        VkPipelineRasterizationStateCreateInfo rasterizer{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
        rasterizer.lineWidth = 1.0f;
        rasterizer.cullMode = VK_CULL_MODE_NONE;
        rasterizer.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;

        VkPipelineMultisampleStateCreateInfo multisampling{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
        multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

        VkPipelineDepthStencilStateCreateInfo depthStencil{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
        depthStencil.depthTestEnable = VK_FALSE;
        depthStencil.depthWriteEnable = VK_FALSE;

        VkPipelineColorBlendAttachmentState colorBlendAttachment{};
        colorBlendAttachment.colorWriteMask = 0xF;
        colorBlendAttachment.blendEnable = VK_FALSE;
        VkPipelineColorBlendStateCreateInfo colorBlending{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
        colorBlending.attachmentCount = 1;
        colorBlending.pAttachments = &colorBlendAttachment;

        VkDynamicState dynamicStates[] = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };
        VkPipelineDynamicStateCreateInfo dynamicState{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
        dynamicState.dynamicStateCount = 2;
        dynamicState.pDynamicStates = dynamicStates;

        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_VERTEX_BIT;
        pushRange.size = sizeof(Matrix4x4);

        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &sceneDepthSetLayout;
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        VkPipelineRenderingCreateInfo renderingInfo{VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO};
        renderingInfo.colorAttachmentCount = 1;
        renderingInfo.pColorAttachmentFormats = &colorFormat;
        renderingInfo.depthAttachmentFormat = VK_FORMAT_D32_SFLOAT;

        VkGraphicsPipelineCreateInfo pipelineInfo{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        pipelineInfo.pNext = &renderingInfo;
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = shaderStages;
        pipelineInfo.pVertexInputState = &vertexInputInfo;
        pipelineInfo.pInputAssemblyState = &inputAssembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &rasterizer;
        pipelineInfo.pMultisampleState = &multisampling;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &colorBlending;
        pipelineInfo.pDynamicState = &dynamicState;
        pipelineInfo.layout = m_pipelineLayout;

        if (vkCreateGraphicsPipelines(nativeDevice, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &m_pipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create debug line pipeline!");
        }

        vkDestroyShaderModule(nativeDevice, vertShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragShader, nullptr);

        m_lineBuffer = std::make_unique<VulkanBuffer>(
            device,
            m_maxLineBufferSize,
            VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
            VMA_MEMORY_USAGE_CPU_TO_GPU
        );
    }

    void DebugLinePass::Shutdown(VkDevice device) {
        m_lineBuffer.reset();
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
    }

    void DebugLinePass::Render(const RenderFrameContext& ctx, const std::vector<DebugLineVertex>& lines) {
        if (lines.empty() || m_pipeline == VK_NULL_HANDLE || !m_lineBuffer) {
            return;
        }

        size_t bufferSize = lines.size() * sizeof(DebugLineVertex);
        if (bufferSize > m_maxLineBufferSize) {
            bufferSize = m_maxLineBufferSize;
        }

        m_lineBuffer->WriteData(lines.data(), bufferSize);

        vkCmdBindPipeline(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
        vkCmdBindDescriptorSets(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &ctx.sceneDepthSet, 0, nullptr);
        Matrix4x4 lineVP = ctx.viewProj.transpose();
        vkCmdPushConstants(ctx.cmd, m_pipelineLayout, VK_SHADER_STAGE_VERTEX_BIT, 0, sizeof(Matrix4x4), &lineVP);

        VkDeviceSize offsets[] = {0};
        VkBuffer rawBuffer = m_lineBuffer->GetNativeBuffer();
        vkCmdBindVertexBuffers(ctx.cmd, 0, 1, &rawBuffer, offsets);

        uint32_t drawCount = static_cast<uint32_t>(bufferSize / sizeof(DebugLineVertex));
        vkCmdDraw(ctx.cmd, drawCount, 1, 0, 0);
    }

    void EditorGridPass::Initialize(VulkanDevice* device, VkFormat colorFormat, VkDescriptorSetLayout sceneDepthSetLayout) {
        VkDevice nativeDevice = device->GetNativeDevice();
        const std::string shaderDir = LIZERAL_SHADER_DIR;

        VkShaderModule vertShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "editor_grid_vert.spv"));
        VkShaderModule fragShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "editor_grid_frag.spv"));

        VkPipelineShaderStageCreateInfo shaderStages[] = {
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_VERTEX_BIT, vertShader, "main", nullptr},
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_FRAGMENT_BIT, fragShader, "main", nullptr}
        };

        VkPipelineVertexInputStateCreateInfo vertexInputInfo{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
        VkPipelineInputAssemblyStateCreateInfo inputAssembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
        inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

        VkPipelineViewportStateCreateInfo viewportState{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
        viewportState.viewportCount = 1;
        viewportState.scissorCount = 1;

        VkPipelineRasterizationStateCreateInfo rasterizer{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
        rasterizer.lineWidth = 1.0f;
        rasterizer.cullMode = VK_CULL_MODE_NONE;
        rasterizer.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;

        VkPipelineMultisampleStateCreateInfo multisampling{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
        multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

        VkPipelineDepthStencilStateCreateInfo depthStencil{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
        depthStencil.depthTestEnable = VK_FALSE;
        depthStencil.depthWriteEnable = VK_FALSE;

        VkPipelineColorBlendAttachmentState colorBlendAttachment{};
        colorBlendAttachment.colorWriteMask = 0xF;
        colorBlendAttachment.blendEnable = VK_TRUE;
        colorBlendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
        colorBlendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        colorBlendAttachment.colorBlendOp = VK_BLEND_OP_ADD;
        colorBlendAttachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
        colorBlendAttachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        colorBlendAttachment.alphaBlendOp = VK_BLEND_OP_ADD;

        VkPipelineColorBlendStateCreateInfo colorBlending{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
        colorBlending.attachmentCount = 1;
        colorBlending.pAttachments = &colorBlendAttachment;

        VkDynamicState dynamicStates[] = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };
        VkPipelineDynamicStateCreateInfo dynamicState{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
        dynamicState.dynamicStateCount = 2;
        dynamicState.pDynamicStates = dynamicStates;

        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
        pushRange.size = sizeof(EditorGridPushConstants);

        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &sceneDepthSetLayout;
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        VkPipelineRenderingCreateInfo renderingInfo{VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO};
        renderingInfo.colorAttachmentCount = 1;
        renderingInfo.pColorAttachmentFormats = &colorFormat;
        renderingInfo.depthAttachmentFormat = VK_FORMAT_D32_SFLOAT;

        VkGraphicsPipelineCreateInfo pipelineInfo{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        pipelineInfo.pNext = &renderingInfo;
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = shaderStages;
        pipelineInfo.pVertexInputState = &vertexInputInfo;
        pipelineInfo.pInputAssemblyState = &inputAssembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &rasterizer;
        pipelineInfo.pMultisampleState = &multisampling;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &colorBlending;
        pipelineInfo.pDynamicState = &dynamicState;
        pipelineInfo.layout = m_pipelineLayout;

        if (vkCreateGraphicsPipelines(nativeDevice, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &m_pipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create editor grid pipeline!");
        }

        vkDestroyShaderModule(nativeDevice, vertShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragShader, nullptr);
    }

    void EditorGridPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
    }

    void EditorGridPass::Render(const RenderFrameContext& ctx, const EditorViewportOverlayData& overlay) {
        if (!overlay.enabled || !overlay.grid.enabled || m_pipeline == VK_NULL_HANDLE) {
            return;
        }

        const float planeHeight = overlay.grid.planeHeight;
        const float cameraPlaneDistance = std::max(std::fabs(ctx.cameraPos.y - planeHeight), 1.0f);
        const float focusLevel = std::max(std::log(std::max(cameraPlaneDistance * 0.25f, 1.0f)) / std::log(5.0f), 0.0f);
        const float baseSpacing = std::max(overlay.grid.minorSpacing, 0.001f);
        const float fadeDistance = std::max(overlay.grid.fadeDistance, 1.0f) * std::max(std::pow(5.0f, std::floor(focusLevel)), 1.0f);

        vkCmdBindPipeline(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
        vkCmdBindDescriptorSets(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &ctx.sceneDepthSet, 0, nullptr);

        for (int level = 5; level >= 0; --level) {
            const float spacing = baseSpacing * std::pow(5.0f, static_cast<float>(level));
            const float focusDistance = std::fabs(static_cast<float>(level) - focusLevel);
            const float focusWeight = 1.0f - SmoothStep(0.0f, 1.75f, focusDistance);
            const float persistentWeight = 1.0f / (1.0f + focusDistance * focusDistance * 1.35f);
            const float opacity = overlay.grid.minorOpacity * (0.18f * persistentWeight + 1.35f * focusWeight);
            if (opacity <= 0.006f) {
                continue;
            }

            EditorGridPushConstants gridPc{};
            gridPc.viewProj = ctx.viewProj.transpose();
            gridPc.invViewProj = ctx.viewProj.inverse().transpose();
            gridPc.cameraPosAndPlaneHeight = Vector4(ctx.cameraPos, planeHeight);
            gridPc.viewportSizeAndSpacing =
                Vector4(static_cast<float>(ctx.width), static_cast<float>(ctx.height), spacing, spacing * 5.0f);
            gridPc.fadeAndOpacity =
                Vector4(fadeDistance,
                        0.0f,
                        std::clamp(opacity, 0.0f, 0.62f),
                        0.0f);

            vkCmdPushConstants(ctx.cmd, m_pipelineLayout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(EditorGridPushConstants), &gridPc);
            vkCmdDraw(ctx.cmd, 3, 1, 0, 0);
        }
    }

    void EditorAxisPass::Initialize(VulkanDevice* device, VkFormat colorFormat, VkDescriptorSetLayout sceneDepthSetLayout) {
        VkDevice nativeDevice = device->GetNativeDevice();
        const std::string shaderDir = LIZERAL_SHADER_DIR;

        VkShaderModule vertShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "editor_axis_vert.spv"));
        VkShaderModule fragShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "editor_axis_frag.spv"));

        VkPipelineShaderStageCreateInfo shaderStages[] = {
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_VERTEX_BIT, vertShader, "main", nullptr},
            {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_FRAGMENT_BIT, fragShader, "main", nullptr}
        };

        VkPipelineVertexInputStateCreateInfo vertexInputInfo{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
        VkPipelineInputAssemblyStateCreateInfo inputAssembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
        inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

        VkPipelineViewportStateCreateInfo viewportState{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
        viewportState.viewportCount = 1;
        viewportState.scissorCount = 1;

        VkPipelineRasterizationStateCreateInfo rasterizer{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
        rasterizer.lineWidth = 1.0f;
        rasterizer.cullMode = VK_CULL_MODE_NONE;
        rasterizer.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;

        VkPipelineMultisampleStateCreateInfo multisampling{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
        multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

        VkPipelineDepthStencilStateCreateInfo depthStencil{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
        depthStencil.depthTestEnable = VK_TRUE;
        depthStencil.depthWriteEnable = VK_TRUE;
        depthStencil.depthCompareOp = VK_COMPARE_OP_GREATER_OR_EQUAL;

        VkPipelineColorBlendAttachmentState colorBlendAttachment{};
        colorBlendAttachment.colorWriteMask = 0xF;
        colorBlendAttachment.blendEnable = VK_TRUE;
        colorBlendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
        colorBlendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        colorBlendAttachment.colorBlendOp = VK_BLEND_OP_ADD;
        colorBlendAttachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
        colorBlendAttachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        colorBlendAttachment.alphaBlendOp = VK_BLEND_OP_ADD;

        VkPipelineColorBlendStateCreateInfo colorBlending{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
        colorBlending.attachmentCount = 1;
        colorBlending.pAttachments = &colorBlendAttachment;

        VkDynamicState dynamicStates[] = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };
        VkPipelineDynamicStateCreateInfo dynamicState{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
        dynamicState.dynamicStateCount = 2;
        dynamicState.pDynamicStates = dynamicStates;

        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
        pushRange.size = sizeof(EditorAxisPushConstants);

        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &sceneDepthSetLayout;
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        VkPipelineRenderingCreateInfo renderingInfo{VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO};
        renderingInfo.colorAttachmentCount = 1;
        renderingInfo.pColorAttachmentFormats = &colorFormat;
        renderingInfo.depthAttachmentFormat = VK_FORMAT_D32_SFLOAT;

        VkGraphicsPipelineCreateInfo pipelineInfo{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        pipelineInfo.pNext = &renderingInfo;
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = shaderStages;
        pipelineInfo.pVertexInputState = &vertexInputInfo;
        pipelineInfo.pInputAssemblyState = &inputAssembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &rasterizer;
        pipelineInfo.pMultisampleState = &multisampling;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &colorBlending;
        pipelineInfo.pDynamicState = &dynamicState;
        pipelineInfo.layout = m_pipelineLayout;

        if (vkCreateGraphicsPipelines(nativeDevice, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &m_pipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create editor axis pipeline!");
        }

        vkDestroyShaderModule(nativeDevice, vertShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragShader, nullptr);
    }

    void EditorAxisPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
    }

    void EditorAxisPass::Render(const RenderFrameContext& ctx, const EditorViewportOverlayData& overlay) {
        if (!overlay.enabled || !overlay.axes.showAxes || m_pipeline == VK_NULL_HANDLE) {
            return;
        }

        struct ViewAxis {
            Vector3 worldAxis;
            Vector4 color;
            char label;
            float sortDepth;
        };

        std::array<ViewAxis, 3> axes = {{
            { Vector3(1.0f, 0.0f, 0.0f), Vector4(0.92f, 0.06f, 0.05f, 1.0f), 'X', 0.0f },
            { Vector3(0.0f, 1.0f, 0.0f), Vector4(0.08f, 0.78f, 0.12f, 1.0f), 'Y', 0.0f },
            { Vector3(0.0f, 0.0f, -1.0f), Vector4(0.10f, 0.34f, 1.0f, 1.0f), 'Z', 0.0f }
        }};

        const float viewportW = static_cast<float>(std::max(ctx.width, 1u));
        const float viewportH = static_cast<float>(std::max(ctx.height, 1u));
        const float widgetSize = std::clamp(std::min(viewportW, viewportH) * 0.13f, 76.0f, 112.0f);
        const float widgetMargin = 14.0f;
        const float widgetMinX = widgetMargin;
        const float widgetMinY = std::max(widgetMargin, viewportH - widgetMargin - widgetSize);
        const float axisThickness = 1.45f;
        const float arrowThickness = 1.25f;
        const float labelThickness = 1.05f;

        auto toClip = [&](float px, float py) {
            return Vector4((px / viewportW) * 2.0f - 1.0f,
                           (py / viewportH) * 2.0f - 1.0f,
                           1.0f,
                           1.0f);
        };

        Matrix4x4 axisView = ctx.view;
        axisView[0][3] = 0.0f;
        axisView[1][3] = 0.0f;
        axisView[2][3] = -3.15f;
        axisView[3][0] = 0.0f;
        axisView[3][1] = 0.0f;
        axisView[3][2] = 0.0f;
        axisView[3][3] = 1.0f;

        constexpr float pi = 3.14159265358979323846f;
        const float fov = 34.0f * pi / 180.0f;
        const float tanHalfFov = std::tan(fov * 0.5f);
        Matrix4x4 axisProj = Matrix4x4::ZERO;
        axisProj[0][0] = 1.0f / tanHalfFov;
        axisProj[1][1] = -1.0f / tanHalfFov;
        axisProj[2][2] = 0.0f;
        axisProj[2][3] = 0.01f;
        axisProj[3][2] = -1.0f;
        axisProj[3][3] = 0.0f;

        const Matrix4x4 axisViewProj = axisProj * axisView;

        struct ProjectedAxisPoint {
            Vector4 clip;
            float pixelX;
            float pixelY;
            float cameraZ;
            bool valid;
        };

        auto projectAxisPoint = [&](const Vector3& position) {
            ProjectedAxisPoint result{};
            Vector4 cameraSpace = axisView * Vector4(position, 1.0f);
            Vector4 clip = axisViewProj * Vector4(position, 1.0f);
            float w = clip.w;
            if (std::fabs(w) < 1e-5f) {
                w = (w < 0.0f) ? -1e-5f : 1e-5f;
            }

            const float ndcX = clip.x / w;
            const float ndcY = clip.y / w;
            result.pixelX = widgetMinX + (ndcX * 0.5f + 0.5f) * widgetSize;
            result.pixelY = widgetMinY + (ndcY * 0.5f + 0.5f) * widgetSize;
            result.clip = toClip(result.pixelX, result.pixelY);
            result.cameraZ = cameraSpace.z;
            result.valid = clip.w > 0.0f;
            return result;
        };

        const float axisLength = 1.0f;
        for (auto& axis : axes) {
            axis.sortDepth = projectAxisPoint(axis.worldAxis * axisLength).cameraZ;
        }

        std::sort(axes.begin(), axes.end(), [](const ViewAxis& lhs, const ViewAxis& rhs) {
            return lhs.sortDepth < rhs.sortDepth;
        });

        EditorAxisPushConstants axisPc{};
        axisPc.viewProj = Matrix4x4::IDENTITY;
        axisPc.viewportSizeAndThickness = Vector4(viewportW, viewportH, axisThickness, 0.0f);

        vkCmdBindPipeline(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
        vkCmdBindDescriptorSets(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &ctx.sceneDepthSet, 0, nullptr);

        auto drawClipLine = [&](const Vector4& a, const Vector4& b, const Vector4& color, float thickness) {
            axisPc.viewportSizeAndThickness.z = thickness;
            axisPc.axisStart = a;
            axisPc.axisEnd = b;
            axisPc.axisColor = color;
            vkCmdPushConstants(ctx.cmd, m_pipelineLayout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(EditorAxisPushConstants), &axisPc);
            vkCmdDraw(ctx.cmd, 3, 1, 0, 0);
        };

        auto drawScreenLine = [&](float x0, float y0, float x1, float y1, const Vector4& color, float thickness) {
            drawClipLine(toClip(x0, y0), toClip(x1, y1), color, thickness);
        };

        auto drawAxisLine3D = [&](const Vector3& a, const Vector3& b, const Vector4& color, float thickness) {
            ProjectedAxisPoint pa = projectAxisPoint(a);
            ProjectedAxisPoint pb = projectAxisPoint(b);
            if (!pa.valid && !pb.valid) {
                return;
            }
            drawClipLine(pa.clip, pb.clip, color, thickness);
        };

        auto drawLabel = [&](char label, float cx, float cy, const Vector4& color) {
            const float s = 3.6f;
            if (label == 'X') {
                drawScreenLine(cx - s, cy - s, cx + s, cy + s, color, labelThickness);
                drawScreenLine(cx - s, cy + s, cx + s, cy - s, color, labelThickness);
            } else if (label == 'Y') {
                drawScreenLine(cx - s, cy - s, cx, cy, color, labelThickness);
                drawScreenLine(cx + s, cy - s, cx, cy, color, labelThickness);
                drawScreenLine(cx, cy, cx, cy + s, color, labelThickness);
            } else if (label == 'Z') {
                drawScreenLine(cx - s, cy - s, cx + s, cy - s, color, labelThickness);
                drawScreenLine(cx + s, cy - s, cx - s, cy + s, color, labelThickness);
                drawScreenLine(cx - s, cy + s, cx + s, cy + s, color, labelThickness);
            }
        };

        const ProjectedAxisPoint projectedOrigin = projectAxisPoint(Vector3::ZERO);

        for (const auto& axis : axes) {
            const Vector3 axisDir = axis.worldAxis.normalisedCopy();
            const Vector3 end = axisDir * axisLength;
            drawAxisLine3D(Vector3::ZERO, end, axis.color, axisThickness);

            Vector3 arrowSide = axisDir.crossProduct(ctx.cameraForward);
            if (arrowSide.squaredLength() < 1e-5f) {
                arrowSide = axisDir.crossProduct(ctx.cameraRight);
            }
            if (arrowSide.squaredLength() < 1e-5f) {
                arrowSide = axisDir.crossProduct(ctx.cameraUp);
            }
            arrowSide.normalise();

            const Vector3 arrowBase = end - axisDir * 0.18f;
            drawAxisLine3D(end, arrowBase + arrowSide * 0.075f, axis.color, arrowThickness);
            drawAxisLine3D(end, arrowBase - arrowSide * 0.075f, axis.color, arrowThickness);

            const ProjectedAxisPoint labelPoint = projectAxisPoint(axisDir * 1.15f);
            const float labelDirX = labelPoint.pixelX - projectedOrigin.pixelX;
            const float labelDirY = labelPoint.pixelY - projectedOrigin.pixelY;
            const float labelDirLen = std::sqrt(labelDirX * labelDirX + labelDirY * labelDirY);
            if (labelDirLen > 1e-3f) {
                drawLabel(axis.label,
                          labelPoint.pixelX + labelDirX / labelDirLen * 3.0f,
                          labelPoint.pixelY + labelDirY / labelDirLen * 3.0f,
                          axis.color);
            }
        }

        drawScreenLine(projectedOrigin.pixelX - 2.7f, projectedOrigin.pixelY,
                       projectedOrigin.pixelX + 2.7f, projectedOrigin.pixelY,
                       Vector4(0.94f, 0.94f, 0.90f, 1.0f), 1.35f);
        drawScreenLine(projectedOrigin.pixelX, projectedOrigin.pixelY - 2.7f,
                       projectedOrigin.pixelX, projectedOrigin.pixelY + 2.7f,
                       Vector4(0.94f, 0.94f, 0.90f, 1.0f), 1.35f);
    }

} // namespace Lizeral
