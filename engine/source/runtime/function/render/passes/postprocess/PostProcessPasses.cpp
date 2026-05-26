#include "runtime/function/render/passes/postprocess/PostProcessPasses.h"

#include "runtime/function/render/passes/RenderPassUtils.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/rhi/vulkan/VulkanPipelineBuilder.h"
#include "runtime/function/render/rhi/vulkan/VulkanShaderUtils.h"

#include <array>
#include <vector>

namespace Lizeral {
namespace {

    VkPushConstantRange CreateFullscreenPushRange() {
        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
        pushRange.size = sizeof(FullscreenPassPushConstants);
        return pushRange;
    }

    void DrawFullscreen(
        const RenderFrameContext& ctx,
        VkPipeline pipeline,
        VkPipelineLayout pipelineLayout,
        VkDescriptorSet primarySet,
        const FullscreenPassPushConstants& pushConstants,
        GBufferAttachment* outputs,
        uint32_t outputCount,
        VkDescriptorSet extraSet = VK_NULL_HANDLE
    ) {
        VkRenderingAttachmentInfo attachments[2] = {};
        for (uint32_t i = 0; i < outputCount; ++i) {
            TransitionImageLayout(ctx.cmd, outputs[i].image, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
            attachments[i].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
            attachments[i].imageView = outputs[i].view;
            attachments[i].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
            attachments[i].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
            attachments[i].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        }

        VkRenderingInfo renderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO};
        renderInfo.renderArea = ctx.scissor;
        renderInfo.layerCount = 1;
        renderInfo.colorAttachmentCount = outputCount;
        renderInfo.pColorAttachments = attachments;

        vkCmdBeginRendering(ctx.cmd, &renderInfo);
        vkCmdSetViewport(ctx.cmd, 0, 1, &ctx.viewport);
        vkCmdSetScissor(ctx.cmd, 0, 1, &ctx.scissor);
        vkCmdBindPipeline(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline);
        if (extraSet != VK_NULL_HANDLE) {
            VkDescriptorSet sets[2] = { primarySet, extraSet };
            vkCmdBindDescriptorSets(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout, 0, 2, sets, 0, nullptr);
        } else {
            vkCmdBindDescriptorSets(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout, 0, 1, &primarySet, 0, nullptr);
        }
        vkCmdPushConstants(ctx.cmd, pipelineLayout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(FullscreenPassPushConstants), &pushConstants);
        vkCmdDraw(ctx.cmd, 3, 1, 0, 0);
        vkCmdEndRendering(ctx.cmd);

        for (uint32_t i = 0; i < outputCount; ++i) {
            TransitionImageLayout(ctx.cmd, outputs[i].image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        }
    }

    VkShaderModule LoadFullscreenVertexShader(VkDevice device) {
        return CreateShaderModule(device, ReadShaderFile(std::string(LIZERAL_SHADER_DIR) + "lighting_vert.spv"));
    }

    VkPipeline BuildSingleOutputPipeline(
        VulkanDevice* device,
        VkPipelineLayout layout,
        VkShaderModule vertexShader,
        VkShaderModule fragmentShader,
        VkFormat colorFormat,
        const VkSpecializationInfo* specializationInfo = nullptr
    ) {
        return VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, vertexShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragmentShader, specializationInfo)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false)
            .SetPipelineLayout(layout)
            .Build(device, { colorFormat }, VK_FORMAT_UNDEFINED);
    }

} // namespace

    void LightingPass::Initialize(VulkanDevice* device, VkDescriptorSetLayout bindlessSetLayout, const LightingPassSettings& settings) {
        VkDevice nativeDevice = device->GetNativeDevice();

        for (uint32_t ping = 0; ping < 2; ++ping) {
            VkDescriptorSetLayoutBinding bindings[5] = {};
            for (uint32_t binding = 0; binding < 4; ++binding) {
                bindings[binding].binding = binding;
                bindings[binding].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                bindings[binding].descriptorCount = 1;
                bindings[binding].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }
            bindings[4].binding = 4;
            bindings[4].descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
            bindings[4].descriptorCount = 1;
            bindings[4].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;

            VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            layoutInfo.bindingCount = 5;
            layoutInfo.pBindings = bindings;
            vkCreateDescriptorSetLayout(nativeDevice, &layoutInfo, nullptr, &m_layouts[ping]);

            VkDescriptorPoolSize poolSizes[2] = {
                { VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 4 },
                { VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR, 1 }
            };
            VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            poolInfo.poolSizeCount = 2;
            poolInfo.pPoolSizes = poolSizes;
            poolInfo.maxSets = 1;
            vkCreateDescriptorPool(nativeDevice, &poolInfo, nullptr, &m_pools[ping]);

            VkDescriptorSetAllocateInfo allocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            allocInfo.descriptorPool = m_pools[ping];
            allocInfo.descriptorSetCount = 1;
            allocInfo.pSetLayouts = &m_layouts[ping];
            vkAllocateDescriptorSets(nativeDevice, &allocInfo, &m_sets[ping]);
        }

        VkPushConstantRange pushRange = CreateFullscreenPushRange();
        VkDescriptorSetLayout setLayouts[2] = { m_layouts[0], bindlessSetLayout };
        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 2;
        pipelineLayoutInfo.pSetLayouts = setLayouts;
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        struct LightingSpecializationData {
            int32_t giQualityLevel;
            int32_t shadowQualityLevel;
        } specializationData{ settings.giQualityLevel, settings.shadowQualityLevel };

        VkSpecializationMapEntry specEntries[2] = {};
        specEntries[0].constantID = 0;
        specEntries[0].offset = offsetof(LightingSpecializationData, giQualityLevel);
        specEntries[0].size = sizeof(int32_t);
        specEntries[1].constantID = 1;
        specEntries[1].offset = offsetof(LightingSpecializationData, shadowQualityLevel);
        specEntries[1].size = sizeof(int32_t);

        VkSpecializationInfo specInfo{};
        specInfo.mapEntryCount = 2;
        specInfo.pMapEntries = specEntries;
        specInfo.dataSize = sizeof(LightingSpecializationData);
        specInfo.pData = &specializationData;

        const std::string shaderDir = LIZERAL_SHADER_DIR;
        VkShaderModule vertexShader = LoadFullscreenVertexShader(nativeDevice);
        VkShaderModule fragmentShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "lighting_uber.spv"));

        m_pipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, vertexShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragmentShader, &specInfo)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false)
            .AddColorBlendAttachment(false)
            .SetPipelineLayout(m_pipelineLayout)
            .Build(device, { VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16B16A16_SFLOAT }, VK_FORMAT_UNDEFINED);

        vkDestroyShaderModule(nativeDevice, vertexShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragmentShader, nullptr);
    }

    void LightingPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
        for (uint32_t i = 0; i < 2; ++i) {
            if (m_pools[i]) {
                vkDestroyDescriptorPool(device, m_pools[i], nullptr);
                m_pools[i] = VK_NULL_HANDLE;
            }
            if (m_layouts[i]) {
                vkDestroyDescriptorSetLayout(device, m_layouts[i], nullptr);
                m_layouts[i] = VK_NULL_HANDLE;
            }
            m_sets[i] = VK_NULL_HANDLE;
        }
    }

    void LightingPass::UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot, VkAccelerationStructureKHR tlasHandle) {
        VkDescriptorImageInfo gInfos[4] = {
            { targets.GetSampler(), targets.albedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.normalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.depth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.velocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
        };

        VkWriteDescriptorSet writes[4] = {};
        for (uint32_t binding = 0; binding < 4; ++binding) {
            writes[binding].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[binding].dstSet = m_sets[frameSlot];
            writes[binding].dstBinding = binding;
            writes[binding].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            writes[binding].descriptorCount = 1;
            writes[binding].pImageInfo = &gInfos[binding];
        }
        vkUpdateDescriptorSets(device->GetNativeDevice(), 4, writes, 0, nullptr);

        UpdateAccelerationStructureDescriptor(device, frameSlot, tlasHandle);
    }

    void LightingPass::UpdateAccelerationStructureDescriptor(VulkanDevice* device, uint32_t frameSlot, VkAccelerationStructureKHR tlasHandle) {
        if (tlasHandle == VK_NULL_HANDLE || m_sets[frameSlot] == VK_NULL_HANDLE) {
            return;
        }

        VkWriteDescriptorSetAccelerationStructureKHR asWrite{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET_ACCELERATION_STRUCTURE_KHR};
        asWrite.accelerationStructureCount = 1;
        asWrite.pAccelerationStructures = &tlasHandle;

        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = m_sets[frameSlot];
        write.dstBinding = 4;
        write.dstArrayElement = 0;
        write.descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
        write.descriptorCount = 1;
        write.pNext = &asWrite;
        vkUpdateDescriptorSets(device->GetNativeDevice(), 1, &write, 0, nullptr);
    }

    void LightingPass::Render(const RenderFrameContext& ctx, RenderFrameTargets& targets, uint32_t frameSlot, VkDescriptorSet bindlessSet, const FullscreenPassPushConstants& pushConstants) {
        GBufferAttachment outputs[2] = { targets.directLight, targets.noisyGI };
        DrawFullscreen(ctx, m_pipeline, m_pipelineLayout, m_sets[frameSlot], pushConstants, outputs, 2, bindlessSet);
    }

    void SvgfTemporalPass::Initialize(VulkanDevice* device) {
        VkDevice nativeDevice = device->GetNativeDevice();
        for (uint32_t ping = 0; ping < 2; ++ping) {
            VkDescriptorSetLayoutBinding bindings[6] = {};
            for (uint32_t binding = 0; binding < 6; ++binding) {
                bindings[binding].binding = binding;
                bindings[binding].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                bindings[binding].descriptorCount = 1;
                bindings[binding].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }
            VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            layoutInfo.bindingCount = 6;
            layoutInfo.pBindings = bindings;
            vkCreateDescriptorSetLayout(nativeDevice, &layoutInfo, nullptr, &m_layouts[ping]);

            VkDescriptorPoolSize poolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 6};
            VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            poolInfo.poolSizeCount = 1;
            poolInfo.pPoolSizes = &poolSize;
            poolInfo.maxSets = 1;
            vkCreateDescriptorPool(nativeDevice, &poolInfo, nullptr, &m_pools[ping]);

            VkDescriptorSetAllocateInfo allocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            allocInfo.descriptorPool = m_pools[ping];
            allocInfo.descriptorSetCount = 1;
            allocInfo.pSetLayouts = &m_layouts[ping];
            vkAllocateDescriptorSets(nativeDevice, &allocInfo, &m_sets[ping]);
        }

        VkPushConstantRange pushRange = CreateFullscreenPushRange();
        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &m_layouts[0];
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        const std::string shaderDir = LIZERAL_SHADER_DIR;
        VkShaderModule vertexShader = LoadFullscreenVertexShader(nativeDevice);
        VkShaderModule fragmentShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "svgf_temporal_frag.spv"));
        m_pipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, vertexShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragmentShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false)
            .AddColorBlendAttachment(false)
            .SetPipelineLayout(m_pipelineLayout)
            .Build(device, { VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16_SFLOAT }, VK_FORMAT_UNDEFINED);
        vkDestroyShaderModule(nativeDevice, vertexShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragmentShader, nullptr);
    }

    void SvgfTemporalPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
        for (uint32_t i = 0; i < 2; ++i) {
            if (m_pools[i]) {
                vkDestroyDescriptorPool(device, m_pools[i], nullptr);
                m_pools[i] = VK_NULL_HANDLE;
            }
            if (m_layouts[i]) {
                vkDestroyDescriptorSetLayout(device, m_layouts[i], nullptr);
                m_layouts[i] = VK_NULL_HANDLE;
            }
            m_sets[i] = VK_NULL_HANDLE;
        }
    }

    void SvgfTemporalPass::UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot) {
        uint32_t pong = (frameSlot + 1) % 2;
        VkDescriptorImageInfo infos[6] = {
            { targets.GetSampler(), targets.noisyGI.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.normalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.depth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.velocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.giHistory[pong].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.momentsHistory[pong].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
        };

        VkWriteDescriptorSet writes[6] = {};
        for (uint32_t binding = 0; binding < 6; ++binding) {
            writes[binding].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[binding].dstSet = m_sets[frameSlot];
            writes[binding].dstBinding = binding;
            writes[binding].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            writes[binding].descriptorCount = 1;
            writes[binding].pImageInfo = &infos[binding];
        }
        vkUpdateDescriptorSets(device->GetNativeDevice(), 6, writes, 0, nullptr);
    }

    void SvgfTemporalPass::Render(const RenderFrameContext& ctx, RenderFrameTargets& targets, uint32_t frameSlot, const FullscreenPassPushConstants& pushConstants) {
        GBufferAttachment outputs[2] = { targets.giHistory[frameSlot], targets.momentsHistory[frameSlot] };
        DrawFullscreen(ctx, m_pipeline, m_pipelineLayout, m_sets[frameSlot], pushConstants, outputs, 2);
    }

    void SvgfATrousPass::Initialize(VulkanDevice* device) {
        VkDevice nativeDevice = device->GetNativeDevice();
        for (uint32_t ping = 0; ping < 2; ++ping) {
            VkDescriptorSetLayoutBinding bindings[4] = {};
            for (uint32_t binding = 0; binding < 4; ++binding) {
                bindings[binding].binding = binding;
                bindings[binding].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                bindings[binding].descriptorCount = 1;
                bindings[binding].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }
            VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            layoutInfo.bindingCount = 4;
            layoutInfo.pBindings = bindings;
            vkCreateDescriptorSetLayout(nativeDevice, &layoutInfo, nullptr, &m_layouts[ping]);

            VkDescriptorPoolSize poolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 16};
            VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            poolInfo.poolSizeCount = 1;
            poolInfo.pPoolSizes = &poolSize;
            poolInfo.maxSets = 4;
            vkCreateDescriptorPool(nativeDevice, &poolInfo, nullptr, &m_pools[ping]);

            std::vector<VkDescriptorSetLayout> layouts(4, m_layouts[ping]);
            VkDescriptorSetAllocateInfo allocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            allocInfo.descriptorPool = m_pools[ping];
            allocInfo.descriptorSetCount = 4;
            allocInfo.pSetLayouts = layouts.data();
            vkAllocateDescriptorSets(nativeDevice, &allocInfo, m_sets[ping]);
        }

        VkPushConstantRange pushRange = CreateFullscreenPushRange();
        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &m_layouts[0];
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        const std::string shaderDir = LIZERAL_SHADER_DIR;
        VkShaderModule vertexShader = LoadFullscreenVertexShader(nativeDevice);
        VkShaderModule fragmentShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "svgf_atrous_frag.spv"));
        m_pipeline = BuildSingleOutputPipeline(device, m_pipelineLayout, vertexShader, fragmentShader, VK_FORMAT_R16G16B16A16_SFLOAT);
        vkDestroyShaderModule(nativeDevice, vertexShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragmentShader, nullptr);
    }

    void SvgfATrousPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
        for (uint32_t i = 0; i < 2; ++i) {
            if (m_pools[i]) {
                vkDestroyDescriptorPool(device, m_pools[i], nullptr);
                m_pools[i] = VK_NULL_HANDLE;
            }
            if (m_layouts[i]) {
                vkDestroyDescriptorSetLayout(device, m_layouts[i], nullptr);
                m_layouts[i] = VK_NULL_HANDLE;
            }
            for (uint32_t iter = 0; iter < 4; ++iter) {
                m_sets[i][iter] = VK_NULL_HANDLE;
            }
        }
    }

    void SvgfATrousPass::UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot) {
        VkImageView inputs[4] = {
            targets.giHistory[frameSlot].view,
            targets.denoisedGITemp.view,
            targets.denoisedGI.view,
            targets.denoisedGITemp.view
        };

        for (uint32_t iter = 0; iter < 4; ++iter) {
            VkDescriptorImageInfo infos[4] = {
                { targets.GetSampler(), inputs[iter], VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { targets.GetSampler(), targets.momentsHistory[frameSlot].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { targets.GetSampler(), targets.normalRoughness.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
                { targets.GetSampler(), targets.depth.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
            };

            VkWriteDescriptorSet writes[4] = {};
            for (uint32_t binding = 0; binding < 4; ++binding) {
                writes[binding].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                writes[binding].dstSet = m_sets[frameSlot][iter];
                writes[binding].dstBinding = binding;
                writes[binding].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                writes[binding].descriptorCount = 1;
                writes[binding].pImageInfo = &infos[binding];
            }
            vkUpdateDescriptorSets(device->GetNativeDevice(), 4, writes, 0, nullptr);
        }
    }

    void SvgfATrousPass::Render(const RenderFrameContext& ctx, RenderFrameTargets& targets, uint32_t frameSlot, FullscreenPassPushConstants pushConstants) {
        int stepSizes[] = { 1, 2, 4, 8 };
        GBufferAttachment* currentOutput = &targets.denoisedGITemp;

        for (uint32_t iter = 0; iter < 4; ++iter) {
            pushConstants.stepSize = stepSizes[iter];
            DrawFullscreen(ctx, m_pipeline, m_pipelineLayout, m_sets[frameSlot][iter], pushConstants, currentOutput, 1);
            currentOutput = (currentOutput == &targets.denoisedGITemp) ? &targets.denoisedGI : &targets.denoisedGITemp;
        }
    }

    void TaaPass::Initialize(VulkanDevice* device) {
        VkDevice nativeDevice = device->GetNativeDevice();
        for (uint32_t ping = 0; ping < 2; ++ping) {
            VkDescriptorSetLayoutBinding bindings[5] = {};
            for (uint32_t binding = 0; binding < 5; ++binding) {
                bindings[binding].binding = binding;
                bindings[binding].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
                bindings[binding].descriptorCount = 1;
                bindings[binding].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
            }
            VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            layoutInfo.bindingCount = 5;
            layoutInfo.pBindings = bindings;
            vkCreateDescriptorSetLayout(nativeDevice, &layoutInfo, nullptr, &m_layouts[ping]);

            VkDescriptorPoolSize poolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 5};
            VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            poolInfo.poolSizeCount = 1;
            poolInfo.pPoolSizes = &poolSize;
            poolInfo.maxSets = 1;
            vkCreateDescriptorPool(nativeDevice, &poolInfo, nullptr, &m_pools[ping]);

            VkDescriptorSetAllocateInfo allocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            allocInfo.descriptorPool = m_pools[ping];
            allocInfo.descriptorSetCount = 1;
            allocInfo.pSetLayouts = &m_layouts[ping];
            vkAllocateDescriptorSets(nativeDevice, &allocInfo, &m_sets[ping]);
        }

        VkPushConstantRange pushRange = CreateFullscreenPushRange();
        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &m_layouts[0];
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        const std::string shaderDir = LIZERAL_SHADER_DIR;
        VkShaderModule vertexShader = LoadFullscreenVertexShader(nativeDevice);
        VkShaderModule fragmentShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "taa_frag.spv"));
        m_pipeline = BuildSingleOutputPipeline(device, m_pipelineLayout, vertexShader, fragmentShader, VK_FORMAT_R16G16B16A16_SFLOAT);
        vkDestroyShaderModule(nativeDevice, vertexShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragmentShader, nullptr);
    }

    void TaaPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
        for (uint32_t i = 0; i < 2; ++i) {
            if (m_pools[i]) {
                vkDestroyDescriptorPool(device, m_pools[i], nullptr);
                m_pools[i] = VK_NULL_HANDLE;
            }
            if (m_layouts[i]) {
                vkDestroyDescriptorSetLayout(device, m_layouts[i], nullptr);
                m_layouts[i] = VK_NULL_HANDLE;
            }
            m_sets[i] = VK_NULL_HANDLE;
        }
    }

    void TaaPass::UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot) {
        uint32_t pong = (frameSlot + 1) % 2;
        VkDescriptorImageInfo infos[5] = {
            { targets.GetSampler(), targets.denoisedGI.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.history[pong].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.velocity.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.albedoMetallic.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL },
            { targets.GetSampler(), targets.directLight.view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL }
        };

        VkWriteDescriptorSet writes[5] = {};
        for (uint32_t binding = 0; binding < 5; ++binding) {
            writes[binding].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[binding].dstSet = m_sets[frameSlot];
            writes[binding].dstBinding = binding;
            writes[binding].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            writes[binding].descriptorCount = 1;
            writes[binding].pImageInfo = &infos[binding];
        }
        vkUpdateDescriptorSets(device->GetNativeDevice(), 5, writes, 0, nullptr);
    }

    void TaaPass::Render(const RenderFrameContext& ctx, RenderFrameTargets& targets, uint32_t frameSlot, const FullscreenPassPushConstants& pushConstants) {
        DrawFullscreen(ctx, m_pipeline, m_pipelineLayout, m_sets[frameSlot], pushConstants, &targets.history[frameSlot], 1);
    }

    void FinalBlitPass::Initialize(VulkanDevice* device, VkFormat colorFormat) {
        VkDevice nativeDevice = device->GetNativeDevice();
        for (uint32_t ping = 0; ping < 2; ++ping) {
            VkDescriptorSetLayoutBinding binding{};
            binding.binding = 0;
            binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            binding.descriptorCount = 1;
            binding.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;

            VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
            layoutInfo.bindingCount = 1;
            layoutInfo.pBindings = &binding;
            vkCreateDescriptorSetLayout(nativeDevice, &layoutInfo, nullptr, &m_layouts[ping]);

            VkDescriptorPoolSize poolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1};
            VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
            poolInfo.poolSizeCount = 1;
            poolInfo.pPoolSizes = &poolSize;
            poolInfo.maxSets = 1;
            vkCreateDescriptorPool(nativeDevice, &poolInfo, nullptr, &m_pools[ping]);

            VkDescriptorSetAllocateInfo allocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            allocInfo.descriptorPool = m_pools[ping];
            allocInfo.descriptorSetCount = 1;
            allocInfo.pSetLayouts = &m_layouts[ping];
            vkAllocateDescriptorSets(nativeDevice, &allocInfo, &m_sets[ping]);
        }

        VkPushConstantRange pushRange = CreateFullscreenPushRange();
        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &m_layouts[0];
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        const std::string shaderDir = LIZERAL_SHADER_DIR;
        VkShaderModule vertexShader = LoadFullscreenVertexShader(nativeDevice);
        VkShaderModule fragmentShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "blit_frag.spv"));
        m_pipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_VERTEX_BIT, vertexShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragmentShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_FRONT_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(false)
            .SetPipelineLayout(m_pipelineLayout)
            .Build(device, { colorFormat }, VK_FORMAT_D32_SFLOAT);
        vkDestroyShaderModule(nativeDevice, vertexShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragmentShader, nullptr);
    }

    void FinalBlitPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
        for (uint32_t i = 0; i < 2; ++i) {
            if (m_pools[i]) {
                vkDestroyDescriptorPool(device, m_pools[i], nullptr);
                m_pools[i] = VK_NULL_HANDLE;
            }
            if (m_layouts[i]) {
                vkDestroyDescriptorSetLayout(device, m_layouts[i], nullptr);
                m_layouts[i] = VK_NULL_HANDLE;
            }
            m_sets[i] = VK_NULL_HANDLE;
        }
    }

    void FinalBlitPass::UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets, uint32_t frameSlot) {
        VkDescriptorImageInfo info{ targets.GetSampler(), targets.history[frameSlot].view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };
        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = m_sets[frameSlot];
        write.dstBinding = 0;
        write.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        write.descriptorCount = 1;
        write.pImageInfo = &info;
        vkUpdateDescriptorSets(device->GetNativeDevice(), 1, &write, 0, nullptr);
    }

    void FinalBlitPass::Render(const RenderFrameContext& ctx, uint32_t frameSlot, const FullscreenPassPushConstants& pushConstants) {
        vkCmdBindPipeline(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
        vkCmdBindDescriptorSets(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &m_sets[frameSlot], 0, nullptr);
        vkCmdPushConstants(ctx.cmd, m_pipelineLayout, VK_SHADER_STAGE_FRAGMENT_BIT, 0, sizeof(FullscreenPassPushConstants), &pushConstants);
        vkCmdDraw(ctx.cmd, 3, 1, 0, 0);
    }

} // namespace Lizeral
