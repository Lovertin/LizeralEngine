#include "RenderDrawPacket.h"

#include "runtime/function/render/resources/RenderResourceCache.h"
#include "runtime/function/render/rhi/vulkan/VulkanBLAS.h"
#include "runtime/function/render/scene/RenderSceneSnapshot.h"
#include "runtime/resource/meshletAssetTypes.h"

#include <algorithm>

namespace Lizeral {

    namespace {

        bool HasTransparentMaterial(const VulkanModelResource& resource) {
            for (const auto& material : resource.materialsCpu) {
                if (material.alphaMode == static_cast<int>(MaterialAlphaMode::Blend)) {
                    return true;
                }
            }
            return false;
        }

        VkTransformMatrixKHR ToVkTransform(const Matrix4x4& modelMatrix) {
            VkTransformMatrixKHR vkTransform{};
            vkTransform.matrix[0][0] = modelMatrix[0][0];
            vkTransform.matrix[0][1] = modelMatrix[0][1];
            vkTransform.matrix[0][2] = modelMatrix[0][2];
            vkTransform.matrix[0][3] = modelMatrix[0][3];
            vkTransform.matrix[1][0] = modelMatrix[1][0];
            vkTransform.matrix[1][1] = modelMatrix[1][1];
            vkTransform.matrix[1][2] = modelMatrix[1][2];
            vkTransform.matrix[1][3] = modelMatrix[1][3];
            vkTransform.matrix[2][0] = modelMatrix[2][0];
            vkTransform.matrix[2][1] = modelMatrix[2][1];
            vkTransform.matrix[2][2] = modelMatrix[2][2];
            vkTransform.matrix[2][3] = modelMatrix[2][3];
            return vkTransform;
        }

    } // namespace

    RenderDrawPacket BuildRenderDrawPacket(const RenderSceneSnapshot& scene, RenderResourceCache& resourceCache) {
        RenderDrawPacket packet;
        packet.meshDraws.reserve(scene.modelDraws.size());
        packet.rtInstanceDescs.reserve(scene.modelDraws.size());
        packet.tlasInstances.reserve(scene.modelDraws.size());

        for (const RenderModelDrawItem& sourceDraw : scene.modelDraws) {
            if (!sourceDraw.visible || sourceDraw.modelAssetPath.empty()) {
                continue;
            }

            VulkanModelResource& resource = resourceCache.GetOrLoadModel(sourceDraw.modelAssetPath);
            if (!resource.IsValid() || !resource.blas) {
                continue;
            }

            const uint64_t materialBuffer = resourceCache.ResolveMaterialBufferAddress(
                sourceDraw.entity,
                sourceDraw.modelAssetPath,
                sourceDraw.resourceRevision,
                sourceDraw.materialOverrides,
                resource
            );

            const uint32_t instanceIndex = static_cast<uint32_t>(packet.meshDraws.size());

            RTInstanceDesc instanceDesc{};
            instanceDesc.vertexBuffer = resource.vAddr;
            instanceDesc.indexBuffer = resource.globalIAddr;
            instanceDesc.materialBuffer = materialBuffer;
            instanceDesc.primitiveMaterialIdBuffer = resource.primMatIdAddr;
            instanceDesc.textureOffset = resource.textureOffset;
            instanceDesc.materialCount = resource.materialCount;
            packet.rtInstanceDescs.push_back(instanceDesc);

            VkAccelerationStructureInstanceKHR instance{};
            instance.transform = ToVkTransform(sourceDraw.modelMatrix);
            instance.instanceCustomIndex = instanceIndex;
            instance.mask = 0xFF;
            instance.flags = VK_GEOMETRY_INSTANCE_TRIANGLE_FACING_CULL_DISABLE_BIT_KHR;
            instance.accelerationStructureReference = resource.blas->GetDeviceAddress();
            packet.tlasInstances.push_back(instance);

            RenderMeshDrawItem draw{};
            draw.entity = sourceDraw.entity;
            draw.resource = &resource;
            draw.modelMatrix = sourceDraw.modelMatrix;
            draw.materialBuffer = materialBuffer;
            draw.instanceIndex = instanceIndex;
            packet.meshDraws.push_back(draw);
            packet.entityToInstanceIndex[static_cast<uint32_t>(sourceDraw.entity)] = instanceIndex;

            if (HasTransparentMaterial(resource)) {
                packet.transparentDraws.push_back({
                    instanceIndex,
                    sourceDraw.worldPosition.squaredDistance(scene.camera.position)
                });
            }
        }

        std::sort(packet.transparentDraws.begin(), packet.transparentDraws.end(), [](const RenderTransparentDrawItem& lhs, const RenderTransparentDrawItem& rhs) {
            return lhs.distanceSq > rhs.distanceSq;
        });

        return packet;
    }

    void WriteRenderDrawPacketInstanceData(
        const RenderDrawPacket& packet,
        bool isFirstFrame,
        std::unordered_map<uint32_t, Matrix4x4>& previousModelMatrices,
        VulkanBuffer& globalInstanceBuffer
    ) {
        GPUInstanceData* mappedInstanceData = static_cast<GPUInstanceData*>(globalInstanceBuffer.Map());
        const uint32_t maxInstanceCount = static_cast<uint32_t>(globalInstanceBuffer.GetSize() / sizeof(GPUInstanceData));

        for (uint32_t drawIndex = 0; drawIndex < packet.meshDraws.size() && drawIndex < maxInstanceCount; ++drawIndex) {
            const RenderMeshDrawItem& draw = packet.meshDraws[drawIndex];
            const uint32_t entityId = static_cast<uint32_t>(draw.entity);
            Matrix4x4 prevModel = isFirstFrame ? draw.modelMatrix : previousModelMatrices[entityId];

            mappedInstanceData[drawIndex].Model = draw.modelMatrix.transpose();
            mappedInstanceData[drawIndex].prevModel = prevModel.transpose();
            previousModelMatrices[entityId] = draw.modelMatrix;
        }

        globalInstanceBuffer.Unmap();
    }

} // namespace Lizeral
