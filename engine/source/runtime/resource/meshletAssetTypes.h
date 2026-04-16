#pragma once

#include <cstdint>

namespace Lizeral {

    enum class MaterialAlphaMode : int32_t {
        Opaque = 0,
        Mask = 1,
        Blend = 2
    };

    struct MeshletVertex {
        float pos[3];
        float normal[3];
        float uv[2];
    };

    struct MeshletBounds {
        float center[3];
        float radius;
    };

    struct MeshletDescriptor {
        uint32_t vertexOffset;
        uint32_t triangleOffset;
        uint32_t vertexCount;
        uint32_t triangleCount;
        uint32_t materialID;
    };

    struct MaterialData {
        float baseColorFactor[4];
        float metallicFactor;
        float roughnessFactor;
        int albedoTex;
        int normalTex;
        int ormTex;
        int emissiveTex;
        int alphaMode;
        float alphaCutoff;
        int pad0;
        int pad1;
    };

} // namespace Lizeral
