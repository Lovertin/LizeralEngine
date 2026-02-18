#pragma once
#include "runtime/function/input/input.h"
#include "runtime/function/ecs/registry.h"

namespace Lizeral{
    class TransformComponent;
    class CameraComponent;
    class CameraControlComponent;
}

namespace Lizeral{
    class CameraControlSystem{

    public:
        CameraControlSystem();
        virtual ~CameraControlSystem();

        void Tick(Registry& registry);

    private:
        void UpdateCameraDir(Input& input,TransformComponent& trans,CameraControlComponent& cameraController);

        void UpdateCameraPosForFree(Input& input,TransformComponent& trans,CameraControlComponent& cameraController);
    };
}