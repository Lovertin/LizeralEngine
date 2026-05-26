#include "RenderResourceCache.h"

#include "runtime/function/ecs/components/Model/VulkanModelComponent.h"
#include "runtime/function/render/MeshletBuilder/MeshletModelBuilder.h"
#include "runtime/function/render/rhi/vulkan/VulkanBLAS.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandPool.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"

#include <iostream>
#include <stdexcept>

namespace Lizeral {

    void RenderResourceCache::Initialize(VulkanDevice* device, VulkanCommandPool* resourceCommandPool) {
        m_device = device;
        m_resourceCommandPool = resourceCommandPool;
    }

    void RenderResourceCache::SetBindlessDescriptorSet(VkDescriptorSet descriptorSet) {
        m_bindlessDescriptorSet = descriptorSet;
    }

    void RenderResourceCache::Clear() {
        m_modelCache.clear();
        m_globalTextures.clear();
        m_globalImageInfos.clear();
        m_globalTexturePathCache.clear();
        m_overrideMaterialBuffers.clear();
        m_bindlessDescriptorSet = VK_NULL_HANDLE;
    }

    void RenderResourceCache::UpdateBindlessTextureDescriptor(uint32_t textureIndex) {
        if (m_device == nullptr || m_bindlessDescriptorSet == VK_NULL_HANDLE || textureIndex >= m_globalImageInfos.size()) {
            return;
        }

        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = m_bindlessDescriptorSet;
        write.dstBinding = 0;
        write.dstArrayElement = textureIndex;
        write.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        write.descriptorCount = 1;
        write.pImageInfo = &m_globalImageInfos[textureIndex];
        vkUpdateDescriptorSets(m_device->GetNativeDevice(), 1, &write, 0, nullptr);
    }

    uint32_t RenderResourceCache::GetOrLoadTextureIndex(const std::string& texturePath) {
        if (texturePath.empty()) {
            return static_cast<uint32_t>(-1);
        }

        auto cached = m_globalTexturePathCache.find(texturePath);
        if (cached != m_globalTexturePathCache.end()) {
            return cached->second;
        }

        if (m_globalTextures.size() >= kMaxBindlessTextures) {
            throw std::runtime_error("Bindless texture array is full.");
        }

        auto texture = std::make_unique<VulkanTexture>(m_device, m_resourceCommandPool, texturePath);

        VkDescriptorImageInfo info{};
        info.imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        info.imageView = texture->GetImageView();
        info.sampler = texture->GetSampler();

        const uint32_t textureIndex = static_cast<uint32_t>(m_globalTextures.size());
        m_globalTextures.push_back(std::move(texture));
        m_globalImageInfos.push_back(info);
        m_globalTexturePathCache[texturePath] = textureIndex;
        UpdateBindlessTextureDescriptor(textureIndex);
        return textureIndex;
    }

    void RenderResourceCache::ApplyMaterialOverride(MaterialData& targetMaterial, const VulkanMaterialSlotOverride& overrideData) {
        const uint32_t overrideMask = overrideData.materialInstance.overrideMask;
        if ((overrideMask & Resource::MaterialOverride_BaseColor) != 0u) {
            for (int channel = 0; channel < 4; ++channel) {
                targetMaterial.baseColorFactor[channel] = overrideData.materialInstance.factors.baseColorFactor[channel];
            }
        }
        if ((overrideMask & Resource::MaterialOverride_Metallic) != 0u) {
            targetMaterial.metallicFactor = overrideData.materialInstance.factors.metallicFactor;
        }
        if ((overrideMask & Resource::MaterialOverride_Roughness) != 0u) {
            targetMaterial.roughnessFactor = overrideData.materialInstance.factors.roughnessFactor;
        }
        if ((overrideMask & Resource::MaterialOverride_AlbedoTex) != 0u) {
            if (overrideData.materialInstance.textures.albedoTex >= 0) {
                targetMaterial.albedoTex = overrideData.materialInstance.textures.albedoTex;
            }
        }
        if ((overrideMask & Resource::MaterialOverride_NormalTex) != 0u) {
            if (overrideData.materialInstance.textures.normalTex >= 0) {
                targetMaterial.normalTex = overrideData.materialInstance.textures.normalTex;
            }
        }
        if ((overrideMask & Resource::MaterialOverride_OrmTex) != 0u) {
            if (overrideData.materialInstance.textures.ormTex >= 0) {
                targetMaterial.ormTex = overrideData.materialInstance.textures.ormTex;
            }
        }
        if ((overrideMask & Resource::MaterialOverride_EmissiveTex) != 0u) {
            if (overrideData.materialInstance.textures.emissiveTex >= 0) {
                targetMaterial.emissiveTex = overrideData.materialInstance.textures.emissiveTex;
            }
        }

        auto applyTexturePath = [this](const std::string& path, int* targetIndex) {
            if (targetIndex == nullptr) {
                return;
            }

            if (path.empty()) {
                return;
            }

            try {
                *targetIndex = static_cast<int>(GetOrLoadTextureIndex(path));
            } catch (const std::exception& exception) {
                std::cerr << "[RenderResourceCache] WARNING: Failed to load override texture '" << path
                          << "'. Reason: " << exception.what() << std::endl;
            }
        };

        const VulkanTextureOverrideSet& textureOverrides = overrideData.textureOverrides;
        if (!textureOverrides.albedoTexturePath.empty()) {
            applyTexturePath(textureOverrides.albedoTexturePath, &targetMaterial.albedoTex);
        }
        if (!textureOverrides.normalTexturePath.empty()) {
            applyTexturePath(textureOverrides.normalTexturePath, &targetMaterial.normalTex);
        }
        if (!textureOverrides.ormTexturePath.empty()) {
            applyTexturePath(textureOverrides.ormTexturePath, &targetMaterial.ormTex);
        }
        if (!textureOverrides.emissiveTexturePath.empty()) {
            applyTexturePath(textureOverrides.emissiveTexturePath, &targetMaterial.emissiveTex);
        }
    }

    uint64_t RenderResourceCache::ResolveMaterialBufferAddress(Entity entity, const VulkanModelComponent& modelComp, const VulkanModelResource& modelResource) {
        return ResolveMaterialBufferAddress(
            entity,
            modelComp.getModelAssetPath(),
            modelComp.GetResourceRevision(),
            modelComp.GetMaterialOverrides(),
            modelResource
        );
    }

    uint64_t RenderResourceCache::ResolveMaterialBufferAddress(
        Entity entity,
        const std::string& modelAssetPath,
        uint32_t resourceRevision,
        const std::vector<VulkanMaterialSlotOverride>& materialOverrides,
        const VulkanModelResource& modelResource
    ) {
        if (materialOverrides.empty() || modelResource.materialsCpu.empty()) {
            return modelResource.matAddr;
        }

        const uint32_t entityId = static_cast<uint32_t>(entity);
        OverrideMaterialBufferCacheEntry& cacheEntry = m_overrideMaterialBuffers[entityId];
        if (cacheEntry.componentRevision == resourceRevision
            && cacheEntry.modelAssetPath == modelAssetPath
            && cacheEntry.materialBufferAddress != 0) {
            return cacheEntry.materialBufferAddress;
        }

        std::vector<MaterialData> resolvedMaterials = modelResource.materialsCpu;
        for (const VulkanMaterialSlotOverride& overrideData : materialOverrides) {
            if (!overrideData.enabled || resolvedMaterials.empty()) {
                continue;
            }

            uint32_t targetMaterialIndex = overrideData.materialSlotIndex;
            if (overrideData.meshAssetIndex < modelResource.meshAssets.size()) {
                targetMaterialIndex = modelResource.meshAssets[overrideData.meshAssetIndex].materialIndex;
            }

            if (targetMaterialIndex >= resolvedMaterials.size()) {
                continue;
            }

            MaterialData targetMaterial = resolvedMaterials[targetMaterialIndex];
            if (overrideData.materialInstance.baseMaterialIndex < resolvedMaterials.size()
                && (overrideData.materialInstance.baseMaterialIndex != 0u || targetMaterialIndex == 0u)) {
                targetMaterial = resolvedMaterials[overrideData.materialInstance.baseMaterialIndex];
            }

            ApplyMaterialOverride(targetMaterial, overrideData);
            resolvedMaterials[targetMaterialIndex] = targetMaterial;
        }

        cacheEntry.materialBuffer = CreateBDABuffer(resolvedMaterials);
        cacheEntry.materialBufferAddress = cacheEntry.materialBuffer ? cacheEntry.materialBuffer->GetDeviceAddress() : modelResource.matAddr;
        cacheEntry.componentRevision = resourceRevision;
        cacheEntry.modelAssetPath = modelAssetPath;
        return cacheEntry.materialBufferAddress;
    }

    VulkanModelResource& RenderResourceCache::GetOrLoadModel(const std::string& path) {
        auto cached = m_modelCache.find(path);
        if (cached != m_modelCache.end()) {
            return cached->second;
        }

        std::cout << "[RenderingSystem] Loading new model to GPU: " << path << std::endl;

        uint32_t currentTexOffset = static_cast<uint32_t>(m_globalTextures.size());
        MeshletModelBuilder builder;
        if (!builder.LoadAndSliceModel(path, currentTexOffset)) {
            throw std::runtime_error("Failed to load model asset: " + path);
        }

        if (m_globalTextures.size() + builder.GetAllTextures().size() > kMaxBindlessTextures) {
            throw std::runtime_error("Bindless texture array is full.");
        }

        for (const auto& texData : builder.GetAllTextures()) {
            auto texture = std::make_unique<VulkanTexture>(m_device, m_resourceCommandPool, texData.data(), static_cast<int>(texData.size()));

            VkDescriptorImageInfo info{};
            info.imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
            info.imageView = texture->GetImageView();
            info.sampler = texture->GetSampler();

            m_globalImageInfos.push_back(info);
            m_globalTextures.push_back(std::move(texture));
        }

        VulkanModelResource res;
        res.totalMeshlets = static_cast<uint32_t>(builder.GetMeshlets().size());
        res.textureOffset = currentTexOffset;
        res.textureCount = static_cast<uint32_t>(builder.GetAllTextures().size());
        res.materialCount = static_cast<uint32_t>(builder.GetMaterials().size());
        res.materialsCpu = builder.GetMaterials();
        res.materialAssets = builder.GetMaterialAssets();
        res.meshAssets = builder.GetMeshAssets();

        res.vertexBuffer = CreateBDABuffer(builder.GetVertices());
        res.meshletBuffer = CreateBDABuffer(builder.GetMeshlets());
        res.indexBuffer = CreateBDABuffer(builder.GetMicroIndices());
        res.boundsBuffer = CreateBDABuffer(builder.GetBounds());
        res.materialBuffer = CreateBDABuffer(builder.GetMaterials());

        res.vAddr = res.vertexBuffer ? res.vertexBuffer->GetDeviceAddress() : 0;
        res.mAddr = res.meshletBuffer ? res.meshletBuffer->GetDeviceAddress() : 0;
        res.iAddr = res.indexBuffer ? res.indexBuffer->GetDeviceAddress() : 0;
        res.bAddr = res.boundsBuffer ? res.boundsBuffer->GetDeviceAddress() : 0;
        res.matAddr = res.materialBuffer ? res.materialBuffer->GetDeviceAddress() : 0;

        const auto& vertices = builder.GetVertices();
        const auto& microIndices = builder.GetMicroIndices();
        const auto& meshlets = builder.GetMeshlets();

        std::vector<uint32_t> globalIndices;
        std::vector<uint32_t> primitiveMaterialIds;
        globalIndices.reserve(microIndices.size());
        primitiveMaterialIds.reserve(microIndices.size() / 3);
        for (const auto& meshlet : meshlets) {
            for (uint32_t i = 0; i < meshlet.triangleCount * 3; i++) {
                globalIndices.push_back(meshlet.vertexOffset + microIndices[meshlet.triangleOffset + i]);
            }
            for (uint32_t tri = 0; tri < meshlet.triangleCount; ++tri) {
                primitiveMaterialIds.push_back(meshlet.materialID);
            }
        }

        uint32_t vertexCount = static_cast<uint32_t>(vertices.size());
        uint32_t indexCount = static_cast<uint32_t>(globalIndices.size());
        uint32_t vertexStride = vertices.empty() ? 0 : static_cast<uint32_t>(sizeof(vertices[0]));

        if (vertexCount > 0 && indexCount > 0) {
            res.globalIndexBuffer = CreateBDABuffer(globalIndices);
            res.globalIAddr = res.globalIndexBuffer->GetDeviceAddress();
            res.primitiveMaterialIdBuffer = CreateBDABuffer(primitiveMaterialIds);
            res.primMatIdAddr = res.primitiveMaterialIdBuffer ? res.primitiveMaterialIdBuffer->GetDeviceAddress() : 0;
            if (!res.primitiveMaterialIdBuffer) {
                res.materialCount = 0;
            }

            std::cout << "[VulkanBLAS] Triggering BLAS build for: " << path << std::endl;
            res.blas = std::make_shared<VulkanBLAS>(
                m_device,
                m_resourceCommandPool,
                res.vAddr,
                vertexCount,
                vertexStride,
                res.globalIAddr,
                indexCount
            );
        }

        if (m_bindlessDescriptorSet != VK_NULL_HANDLE && res.textureCount > 0) {
            VkWriteDescriptorSet write{};
            write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write.dstSet = m_bindlessDescriptorSet;
            write.dstBinding = 0;
            write.dstArrayElement = currentTexOffset;
            write.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            write.descriptorCount = res.textureCount;
            write.pImageInfo = m_globalImageInfos.data() + currentTexOffset;
            vkUpdateDescriptorSets(m_device->GetNativeDevice(), 1, &write, 0, nullptr);
        }

        auto inserted = m_modelCache.emplace(path, std::move(res));
        return inserted.first->second;
    }

} // namespace Lizeral
