#pragma once

#include "runtime/core/math/matrix4.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/vector4.h"
#include "runtime/function/ecs/entity.h"
#include "runtime/function/ecs/components/Model/VulkanModelComponent.h"

#include <string>
#include <vector>

namespace Lizeral {

    class Registry;

    struct RenderCameraSnapshot {
        Matrix4x4 view;
        Matrix4x4 projection;
        Matrix4x4 jitteredProjection;
        Vector3 position;
        Vector3 right {1.0f, 0.0f, 0.0f};
        Vector3 up {0.0f, 1.0f, 0.0f};
        Vector3 forward {0.0f, 0.0f, -1.0f};
        bool valid = false;
    };

    struct RenderDirectionalLightSnapshot {
        Vector3 direction {0.0f, -1.0f, 0.0f};
        Vector3 color {1.0f, 1.0f, 1.0f};
        float intensity = 0.0f;
    };

    struct RenderModelDrawItem {
        Entity entity = null_entity;
        std::string modelAssetPath;
        Matrix4x4 modelMatrix;
        Vector3 worldPosition;
        bool visible = false;
        uint32_t resourceRevision = 0;
        std::vector<VulkanMaterialSlotOverride> materialOverrides;
    };

    struct RenderPointLightSnapshot {
        Vector4 posAndRadius;
        Vector4 colorAndIntensity;
    };

    struct RenderSceneSnapshot {
        RenderCameraSnapshot camera;
        RenderDirectionalLightSnapshot directionalLight;
        std::vector<RenderModelDrawItem> modelDraws;
        std::vector<RenderPointLightSnapshot> pointLights;
    };

    RenderSceneSnapshot ExtractRenderSceneSnapshot(Registry& registry, uint32_t width, uint32_t height, float ndcJitterX, float ndcJitterY);

} // namespace Lizeral
