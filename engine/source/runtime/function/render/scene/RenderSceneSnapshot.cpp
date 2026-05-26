#include "RenderSceneSnapshot.h"

#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Light/DirectionalLightComponent.h"
#include "runtime/function/ecs/components/Light/PointLightComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/registry.h"

namespace Lizeral {

    RenderSceneSnapshot ExtractRenderSceneSnapshot(Registry& registry, uint32_t width, uint32_t height, float ndcJitterX, float ndcJitterY) {
        RenderSceneSnapshot snapshot;

        auto cameraView = registry.view<TransformComponent, CameraComponent>();
        for (auto entity : cameraView) {
            auto& camera = cameraView.get<CameraComponent>(entity);
            auto& transform = cameraView.get<TransformComponent>(entity);

            snapshot.camera.position = transform.getPosition();
            snapshot.camera.right = transform.getRight().normalisedCopy();
            snapshot.camera.up = transform.getUp().normalisedCopy();
            snapshot.camera.forward = transform.getForward().normalisedCopy();

            float aspect = static_cast<float>(width) / static_cast<float>(height > 0 ? height : 1);
            Matrix4x4 baseProj = camera.BuildPerspectiveInfiniteReverseZ(camera.getFov(), aspect, camera.getzNear());
            camera.setProjectionMatrix(baseProj);

            snapshot.camera.view = camera.getViewMatrix();
            snapshot.camera.projection = camera.getProjectionMatrix();
            snapshot.camera.projection[1][1] *= -1.0f;

            snapshot.camera.jitteredProjection = snapshot.camera.projection;
            snapshot.camera.jitteredProjection[0][2] += ndcJitterX;
            snapshot.camera.jitteredProjection[1][2] += ndcJitterY;
            snapshot.camera.valid = true;
            break;
        }

        auto lightView = registry.view<TransformComponent, DirectionLightComponent>();
        for (auto entity : lightView) {
            auto& transform = lightView.get<TransformComponent>(entity);
            auto& light = lightView.get<DirectionLightComponent>(entity);
            snapshot.directionalLight.direction = transform.getForward().normalize();
            snapshot.directionalLight.color = light.getColor();
            snapshot.directionalLight.intensity = light.getIntensity();
            break;
        }

        auto modelView = registry.view<TransformComponent, VulkanModelComponent>();
        for (auto entity : modelView) {
            auto& transform = modelView.get<TransformComponent>(entity);
            auto& modelComp = modelView.get<VulkanModelComponent>(entity);

            if (!modelComp.IsVisible() || !modelComp.HasModelAsset()) {
                continue;
            }

            RenderModelDrawItem item;
            item.entity = entity;
            item.modelAssetPath = modelComp.getModelAssetPath();
            item.modelMatrix = transform.getMatrix();
            item.worldPosition = transform.getPosition();
            item.visible = true;
            item.resourceRevision = modelComp.GetResourceRevision();
            item.materialOverrides = modelComp.GetMaterialOverrides();
            snapshot.modelDraws.push_back(std::move(item));
        }

        auto pointLightView = registry.view<TransformComponent, PointLightComponent>();
        for (auto entity : pointLightView) {
            auto& transform = pointLightView.get<TransformComponent>(entity);
            auto& pointLight = pointLightView.get<PointLightComponent>(entity);

            snapshot.pointLights.push_back({
                Vector4(transform.getPosition().x, transform.getPosition().y, transform.getPosition().z, pointLight.getSpotRadius()),
                Vector4(pointLight.getSpotLightColor().x, pointLight.getSpotLightColor().y, pointLight.getSpotLightColor().z, pointLight.getSpotIntensity())
            });
        }

        return snapshot;
    }

} // namespace Lizeral
