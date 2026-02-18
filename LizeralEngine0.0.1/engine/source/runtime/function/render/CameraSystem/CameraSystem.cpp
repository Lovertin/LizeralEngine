#include "CameraSystem.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/input/input.h"

#include "runtime/core/math/matrix4.h"

namespace Lizeral{
    //脏位恰好是为了服务于system的
    const uint32_t PROJECTION_DIRTY_MASK = 
                CAMERA_DIRTY_FOV | 
                CAMERA_DIRTY_ASPECT | 
                CAMERA_DIRTY_ZNEAR | 
                CAMERA_DIRTY_ZFAR;

    CameraSystem::CameraSystem(){}

    CameraSystem::~CameraSystem(){}

    void CameraSystem::Tick(Registry& registry) {

        // 遍历所有同时拥有 Transform 和 Camera 组件的实体
        auto view = registry.view<TransformComponent, CameraComponent>();

        for (auto entity : view) {
            auto& transform = view.get<TransformComponent>(entity);
            auto& camera = view.get<CameraComponent>(entity);

            // --- 阶段 1: 投影矩阵 (低频) ---
            // 检查你在 Component 里定义的脏标记
            // CAMERA_DIRTY_FOV | CAMERA_DIRTY_ASPECT | CAMERA_DIRTY_ZNEAR | CAMERA_DIRTY_ZFAR
            if (camera.isDirty(PROJECTION_DIRTY_MASK)) { 
                UpdateProjectionMatrix(camera);
                // 注意：在 UpdateProjectionMatrix 内部算完后，记得 camera.cleanDirty();
            }

            // --- 阶段 2: 视图矩阵 (高频) ---
            // 只要 Transform 变了，或者强制每帧更新
            // 优化：其实也可以检查 Transform 的 dirty，但相机通常每帧都在动
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
        // 1. 获取相机的世界变换矩阵 (Position + Rotation)
        Matrix4x4 worldMatrix = transform.getMatrix(); // 假设 Transform 组件有这个 helper

        // 2. 视图矩阵 = 世界矩阵的逆矩阵
        // 相机向右移 = 世界向左移
        Matrix4x4 viewMatrix = worldMatrix.inverse();

        // 3. 写入组件缓存
        camera.setViewMatrix(viewMatrix);
    }
}