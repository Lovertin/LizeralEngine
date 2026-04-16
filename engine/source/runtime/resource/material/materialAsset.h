#pragma once

#include "runtime/resource/meshletAssetTypes.h"

#include <cstdint>
#include <string>

namespace Lizeral::Resource {

    enum class MaterialShadingModel : uint8_t {
        PBRMetallicRoughness = 0,
        Unlit = 1,
        Custom = 255
    };

    struct MaterialTextureSlots {
        int albedoTex {-1};
        int normalTex {-1};
        int ormTex {-1};
        int emissiveTex {-1};
    };

    struct MaterialFactors {
        float baseColorFactor[4] {1.0f, 1.0f, 1.0f, 1.0f};
        float metallicFactor {1.0f};
        float roughnessFactor {1.0f};
        float emissiveFactor[3] {0.0f, 0.0f, 0.0f};
    };

    struct MaterialAsset {
        std::string name;
        std::string shaderTag {"pbr_metallic_roughness"};
        MaterialShadingModel shadingModel {MaterialShadingModel::PBRMetallicRoughness};
        MaterialAlphaMode alphaMode {MaterialAlphaMode::Opaque};
        float alphaCutoff {0.5f};
        MaterialFactors factors {};
        MaterialTextureSlots textures {};
    };

    enum MaterialOverrideBits : uint32_t {
        MaterialOverride_None = 0,
        MaterialOverride_BaseColor = 1u << 0,
        MaterialOverride_Metallic = 1u << 1,
        MaterialOverride_Roughness = 1u << 2,
        MaterialOverride_AlbedoTex = 1u << 3,
        MaterialOverride_NormalTex = 1u << 4,
        MaterialOverride_OrmTex = 1u << 5,
        MaterialOverride_EmissiveTex = 1u << 6
    };

    struct MaterialInstance {
        uint32_t baseMaterialIndex {0};
        uint32_t overrideMask {MaterialOverride_None};
        MaterialFactors factors {};
        MaterialTextureSlots textures {};
    };

    inline MaterialAsset MakeMaterialAssetFromGpuMaterial(const MaterialData& material, const std::string& name = {}) {
        MaterialAsset out;
        out.name = name;
        out.factors.baseColorFactor[0] = material.baseColorFactor[0];
        out.factors.baseColorFactor[1] = material.baseColorFactor[1];
        out.factors.baseColorFactor[2] = material.baseColorFactor[2];
        out.factors.baseColorFactor[3] = material.baseColorFactor[3];
        out.factors.metallicFactor = material.metallicFactor;
        out.factors.roughnessFactor = material.roughnessFactor;
        out.alphaMode = static_cast<MaterialAlphaMode>(material.alphaMode);
        out.alphaCutoff = material.alphaCutoff;

        out.textures.albedoTex = material.albedoTex;
        out.textures.normalTex = material.normalTex;
        out.textures.ormTex = material.ormTex;
        out.textures.emissiveTex = material.emissiveTex;
        return out;
    }

} // namespace Lizeral::Resource
