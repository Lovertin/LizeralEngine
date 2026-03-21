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
        void UpdateProjectionMatrix(CameraComponent& camera);

        void UpdateViewMatrix(const TransformComponent& transform, CameraComponent& camera);
        
    };
}