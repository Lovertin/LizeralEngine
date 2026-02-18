#pragma once
#include "runtime/function/ecs/registry.h"

namespace Lizeral {
    class TransformComponent;
    class CameraComponent;
}

namespace Lizeral{
    class CameraSystem{

    public:
        CameraSystem();
        virtual ~CameraSystem();

        void Tick(Registry& registry);

    private:
        // 如果fov\aspect\zNear\zFar dirty则计算proj
        void UpdateProjectionMatrix(CameraComponent& camera);

        // 每帧根据位置和旋转计算view矩阵
        void UpdateViewMatrix(const TransformComponent& transform, CameraComponent& camera);
        
    };
}