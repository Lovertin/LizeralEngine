#include "runtime/function/render/passes/mesh/MeshRenderPasses.h"

#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystemTypes.h"
#include "runtime/function/render/draw/RenderDrawPacket.h"
#include "runtime/function/render/passes/RenderPassUtils.h"
#include "runtime/function/render/resources/RenderResourceCache.h"
#include "runtime/function/render/rhi/vulkan/VulkanDescriptorBuilder.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/rhi/vulkan/VulkanPipelineBuilder.h"
#include "runtime/function/render/rhi/vulkan/VulkanShaderUtils.h"
#include "runtime/function/render/rhi/vulkan/VulkanTLAS.h"
#include "runtime/function/render/targets/RenderFrameTargets.h"

#include <stdexcept>

namespace Lizeral {
namespace {

    VkPushConstantRange CreateMeshPushRange() {
        VkPushConstantRange pushRange{};
        pushRange.stageFlags = VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT;
        pushRange.size = sizeof(MeshPassPushConstants);
        return pushRange;
    }

    void FillMeshPushConstants(
        MeshPassPushConstants& pushData,
        const RenderMeshDrawItem& draw,
        uint64_t frameDataAddr,
        uint64_t baseInstanceAddr
    ) {
        const VulkanModelResource& res = *draw.resource;
        pushData.frameDataAddr = frameDataAddr;
        pushData.modelDataAddr = baseInstanceAddr + (draw.instanceIndex * sizeof(GPUInstanceData));
        pushData.vertexBuffer = res.vAddr;
        pushData.meshletBuffer = res.mAddr;
        pushData.indexBuffer = res.iAddr;
        pushData.boundsBuffer = res.bAddr;
        pushData.materialBuffer = draw.materialBuffer;
        pushData.totalMeshlets = res.totalMeshlets;
        pushData.textureOffset = 0;
    }

} // namespace

    void TlasBuildPass::Render(const RenderFrameContext& ctx, VulkanTLAS& tlas, uint32_t frameSlot, const RenderDrawPacket& packet) {
        if (packet.tlasInstances.empty()) {
            return;
        }

        tlas.Build(ctx.cmd, frameSlot, packet.tlasInstances, false);

        VkMemoryBarrier memoryBarrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
        memoryBarrier.srcAccessMask = VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR;
        memoryBarrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        vkCmdPipelineBarrier(
            ctx.cmd,
            VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR,
            VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT,
            0,
            1,
            &memoryBarrier,
            0,
            nullptr,
            0,
            nullptr
        );
    }

    void GBufferPass::Initialize(VulkanDevice* device, RenderResourceCache& resourceCache, const RenderFrameTargets& targets) {
        VkDevice nativeDevice = device->GetNativeDevice();
        const std::string shaderDir = LIZERAL_SHADER_DIR;

        VulkanDescriptorBuilder bindlessBuilder;
        VkDescriptorImageInfo dummyInfo = {
            targets.GetSampler(),
            targets.albedoMetallic.view,
            VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL
        };
        const auto& globalImageInfos = resourceCache.GetGlobalImageInfos();
        bindlessBuilder.BindImageArray(
            0,
            globalImageInfos.empty() ? &dummyInfo : globalImageInfos.data(),
            RenderResourceCache::kMaxBindlessTextures,
            globalImageInfos.empty() ? 1 : static_cast<uint32_t>(globalImageInfos.size()),
            VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
            VK_SHADER_STAGE_FRAGMENT_BIT,
            true
        );
        bindlessBuilder.Build(device, m_bindlessSetLayout, m_bindlessPool, m_bindlessSet);
        resourceCache.SetBindlessDescriptorSet(m_bindlessSet);

        VkPushConstantRange pushRange = CreateMeshPushRange();
        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &m_bindlessSetLayout;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        VkShaderModule taskShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "car_task.spv"));
        VkShaderModule meshShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "car_mesh.spv"));
        VkShaderModule fragShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "car_frag.spv"));

        m_pipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_TASK_BIT_EXT, taskShader)
            .AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(true, true, VK_COMPARE_OP_GREATER_OR_EQUAL)
            .DisableColorBlendAttachments(3)
            .SetPipelineLayout(m_pipelineLayout)
            .Build(device, { VK_FORMAT_R8G8B8A8_UNORM, VK_FORMAT_R16G16B16A16_SFLOAT, VK_FORMAT_R16G16_SFLOAT }, VK_FORMAT_D32_SFLOAT);

        vkDestroyShaderModule(nativeDevice, taskShader, nullptr);
        vkDestroyShaderModule(nativeDevice, meshShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragShader, nullptr);
    }

    void GBufferPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
        if (m_bindlessPool) {
            vkDestroyDescriptorPool(device, m_bindlessPool, nullptr);
            m_bindlessPool = VK_NULL_HANDLE;
        }
        if (m_bindlessSetLayout) {
            vkDestroyDescriptorSetLayout(device, m_bindlessSetLayout, nullptr);
            m_bindlessSetLayout = VK_NULL_HANDLE;
        }
        m_bindlessSet = VK_NULL_HANDLE;
    }

    void GBufferPass::Render(
        const RenderFrameContext& ctx,
        RenderFrameTargets& targets,
        const RenderDrawPacket& packet,
        uint64_t frameDataAddr,
        uint64_t baseInstanceAddr,
        bool firstFrame,
        PFN_vkCmdDrawMeshTasksEXT cmdDrawMeshTasks
    ) {
        if (m_pipeline == VK_NULL_HANDLE || cmdDrawMeshTasks == nullptr) {
            return;
        }

        VkImageLayout currentLayout = firstFrame ? VK_IMAGE_LAYOUT_UNDEFINED : VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        TransitionImageLayout(ctx.cmd, targets.albedoMetallic.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        TransitionImageLayout(ctx.cmd, targets.normalRoughness.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        TransitionImageLayout(ctx.cmd, targets.depth.image, currentLayout, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);
        TransitionImageLayout(ctx.cmd, targets.velocity.image, currentLayout, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);

        VkRenderingAttachmentInfo colorAttachments[3] = {};
        colorAttachments[0].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
        colorAttachments[0].imageView = targets.albedoMetallic.view;
        colorAttachments[0].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        colorAttachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        colorAttachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        colorAttachments[0].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};
        colorAttachments[1].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
        colorAttachments[1].imageView = targets.normalRoughness.view;
        colorAttachments[1].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        colorAttachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        colorAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        colorAttachments[1].clearValue.color = {0.0f, 0.0f, 0.0f, 1.0f};
        colorAttachments[2].sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
        colorAttachments[2].imageView = targets.velocity.view;
        colorAttachments[2].imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        colorAttachments[2].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        colorAttachments[2].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        colorAttachments[2].clearValue.color = {0.0f, 0.0f, 0.0f, 0.0f};

        VkRenderingAttachmentInfo depthAttachment{VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO};
        depthAttachment.imageView = targets.depth.view;
        depthAttachment.imageLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
        depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        depthAttachment.clearValue.depthStencil = {0.0f, 0};

        VkRenderingInfo renderInfo{VK_STRUCTURE_TYPE_RENDERING_INFO};
        renderInfo.renderArea = ctx.scissor;
        renderInfo.layerCount = 1;
        renderInfo.colorAttachmentCount = 3;
        renderInfo.pColorAttachments = colorAttachments;
        renderInfo.pDepthAttachment = &depthAttachment;

        vkCmdBeginRendering(ctx.cmd, &renderInfo);
        vkCmdSetViewport(ctx.cmd, 0, 1, &ctx.viewport);
        vkCmdSetScissor(ctx.cmd, 0, 1, &ctx.scissor);
        vkCmdBindPipeline(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
        vkCmdBindDescriptorSets(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &m_bindlessSet, 0, nullptr);

        for (const RenderMeshDrawItem& draw : packet.meshDraws) {
            if (draw.resource == nullptr || !draw.resource->IsValid()) {
                continue;
            }

            MeshPassPushConstants pushData{};
            FillMeshPushConstants(pushData, draw, frameDataAddr, baseInstanceAddr);
            vkCmdPushConstants(
                ctx.cmd,
                m_pipelineLayout,
                VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT,
                0,
                sizeof(MeshPassPushConstants),
                &pushData
            );
            cmdDrawMeshTasks(ctx.cmd, (draw.resource->totalMeshlets + 63) / 64, 1, 1);
        }

        vkCmdEndRendering(ctx.cmd);

        TransitionImageLayout(ctx.cmd, targets.albedoMetallic.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        TransitionImageLayout(ctx.cmd, targets.normalRoughness.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
        TransitionImageLayout(ctx.cmd, targets.depth.image, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_DEPTH_BIT);
        TransitionImageLayout(ctx.cmd, targets.velocity.image, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_IMAGE_ASPECT_COLOR_BIT);
    }

    void TransparentPass::Initialize(VulkanDevice* device, VkFormat colorFormat, VkDescriptorSetLayout bindlessSetLayout) {
        VkDevice nativeDevice = device->GetNativeDevice();
        const std::string shaderDir = LIZERAL_SHADER_DIR;

        VkDescriptorSetLayoutBinding depthBinding{};
        depthBinding.binding = 0;
        depthBinding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        depthBinding.descriptorCount = 1;
        depthBinding.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;

        VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        layoutInfo.bindingCount = 1;
        layoutInfo.pBindings = &depthBinding;
        vkCreateDescriptorSetLayout(nativeDevice, &layoutInfo, nullptr, &m_sceneDepthSetLayout);

        VkDescriptorPoolSize poolSize{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1};
        VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
        poolInfo.poolSizeCount = 1;
        poolInfo.pPoolSizes = &poolSize;
        poolInfo.maxSets = 1;
        vkCreateDescriptorPool(nativeDevice, &poolInfo, nullptr, &m_sceneDepthPool);

        VkDescriptorSetAllocateInfo allocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        allocInfo.descriptorPool = m_sceneDepthPool;
        allocInfo.descriptorSetCount = 1;
        allocInfo.pSetLayouts = &m_sceneDepthSetLayout;
        vkAllocateDescriptorSets(nativeDevice, &allocInfo, &m_sceneDepthSet);

        VkPushConstantRange pushRange = CreateMeshPushRange();
        VkDescriptorSetLayout setLayouts[2] = { m_sceneDepthSetLayout, bindlessSetLayout };
        VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pipelineLayoutInfo.setLayoutCount = 2;
        pipelineLayoutInfo.pSetLayouts = setLayouts;
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushRange;
        vkCreatePipelineLayout(nativeDevice, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);

        VkShaderModule taskShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "car_task.spv"));
        VkShaderModule meshShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "car_mesh.spv"));
        VkShaderModule fragShader = CreateShaderModule(nativeDevice, ReadShaderFile(shaderDir + "transparent_frag.spv"));

        m_pipeline = VulkanPipelineBuilder()
            .AddShaderStage(VK_SHADER_STAGE_TASK_BIT_EXT, taskShader)
            .AddShaderStage(VK_SHADER_STAGE_MESH_BIT_EXT, meshShader)
            .AddShaderStage(VK_SHADER_STAGE_FRAGMENT_BIT, fragShader)
            .SetInputAssembly(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST)
            .SetRasterization(VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, VK_FRONT_FACE_COUNTER_CLOCKWISE)
            .SetMultisampling(VK_SAMPLE_COUNT_1_BIT)
            .SetDepthStencil(false, false, VK_COMPARE_OP_ALWAYS)
            .AddColorBlendAttachment(true)
            .SetPipelineLayout(m_pipelineLayout)
            .Build(device, { colorFormat }, VK_FORMAT_D32_SFLOAT);

        vkDestroyShaderModule(nativeDevice, taskShader, nullptr);
        vkDestroyShaderModule(nativeDevice, meshShader, nullptr);
        vkDestroyShaderModule(nativeDevice, fragShader, nullptr);
    }

    void TransparentPass::Shutdown(VkDevice device) {
        if (m_pipeline) {
            vkDestroyPipeline(device, m_pipeline, nullptr);
            m_pipeline = VK_NULL_HANDLE;
        }
        if (m_pipelineLayout) {
            vkDestroyPipelineLayout(device, m_pipelineLayout, nullptr);
            m_pipelineLayout = VK_NULL_HANDLE;
        }
        if (m_sceneDepthPool) {
            vkDestroyDescriptorPool(device, m_sceneDepthPool, nullptr);
            m_sceneDepthPool = VK_NULL_HANDLE;
        }
        if (m_sceneDepthSetLayout) {
            vkDestroyDescriptorSetLayout(device, m_sceneDepthSetLayout, nullptr);
            m_sceneDepthSetLayout = VK_NULL_HANDLE;
        }
        m_sceneDepthSet = VK_NULL_HANDLE;
    }

    void TransparentPass::UpdateDescriptors(VulkanDevice* device, const RenderFrameTargets& targets) {
        if (m_sceneDepthSet == VK_NULL_HANDLE) {
            return;
        }

        VkDescriptorImageInfo depthInfo {
            targets.GetSampler(),
            targets.depth.view,
            VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL
        };
        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = m_sceneDepthSet;
        write.dstBinding = 0;
        write.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        write.descriptorCount = 1;
        write.pImageInfo = &depthInfo;
        vkUpdateDescriptorSets(device->GetNativeDevice(), 1, &write, 0, nullptr);
    }

    void TransparentPass::Render(
        const RenderFrameContext& ctx,
        const RenderDrawPacket& packet,
        VkDescriptorSet bindlessSet,
        uint64_t frameDataAddr,
        uint64_t baseInstanceAddr,
        PFN_vkCmdDrawMeshTasksEXT cmdDrawMeshTasks
    ) {
        if (packet.transparentDraws.empty() || m_pipeline == VK_NULL_HANDLE || cmdDrawMeshTasks == nullptr) {
            return;
        }

        vkCmdBindPipeline(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
        VkDescriptorSet sets[2] = { m_sceneDepthSet, bindlessSet };
        vkCmdBindDescriptorSets(ctx.cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 2, sets, 0, nullptr);

        for (const RenderTransparentDrawItem& transparentDraw : packet.transparentDraws) {
            if (transparentDraw.meshDrawIndex >= packet.meshDraws.size()) {
                continue;
            }

            const RenderMeshDrawItem& draw = packet.meshDraws[transparentDraw.meshDrawIndex];
            if (draw.resource == nullptr || !draw.resource->IsValid()) {
                continue;
            }

            MeshPassPushConstants pushData{};
            FillMeshPushConstants(pushData, draw, frameDataAddr, baseInstanceAddr);
            vkCmdPushConstants(
                ctx.cmd,
                m_pipelineLayout,
                VK_SHADER_STAGE_TASK_BIT_EXT | VK_SHADER_STAGE_MESH_BIT_EXT | VK_SHADER_STAGE_FRAGMENT_BIT,
                0,
                sizeof(MeshPassPushConstants),
                &pushData
            );
            cmdDrawMeshTasks(ctx.cmd, (draw.resource->totalMeshlets + 63) / 64, 1, 1);
        }
    }

} // namespace Lizeral
