#include "CameraSystem.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/input/input.h"

#include "runtime/core/math/matrix4.h"

namespace Lizeral{
    const uint32_t PROJECTION_DIRTY_MASK = 
                CAMERA_DIRTY_FOV | 
                CAMERA_DIRTY_ASPECT | 
                CAMERA_DIRTY_ZNEAR | 
                CAMERA_DIRTY_ZFAR;

    CameraSystem::CameraSystem(){}

    CameraSystem::~CameraSystem(){}

    void CameraSystem::Tick(Registry& registry) {

        auto view = registry.view<TransformComponent, CameraComponent>();

        for (auto entity : view) {
            auto& transform = view.get<TransformComponent>(entity);
            auto& camera = view.get<CameraComponent>(entity);

            if (camera.isDirty(PROJECTION_DIRTY_MASK)) { 
                UpdateProjectionMatrix(camera);
            }

            UpdateViewMatrix(transform, camera);
        }
    }

    void CameraSystem::UpdateProjectionMatrix(CameraComponent& camera){
        float fov=camera.getFov();
        float aspect=camera.getAspect();
        float zNear=camera.getzNear();
        float zFar=camera.getzFar();

        Matrix4x4 projMat=camera.buildPerspective(fov,aspect,zNear,zFar);
        camera.setProjectionMatrix(projMat);
        camera.clearDirty(PROJECTION_DIRTY_MASK);
    }

    void CameraSystem::UpdateViewMatrix(const TransformComponent& transform, CameraComponent& camera){

        Matrix4x4 worldMatrix = transform.getMatrix(); 

        Matrix4x4 viewMatrix = worldMatrix.inverse();

        camera.setViewMatrix(viewMatrix);
    }
}