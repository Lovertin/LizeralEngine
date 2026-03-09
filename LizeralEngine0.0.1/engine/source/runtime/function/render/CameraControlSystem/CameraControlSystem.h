#pragma once
#include "runtime/function/ecs/registry.h"
#include "runtime/function/input/input.h"
// #include <GLFW/glfw3.h>

namespace Lizeral{
    class TransformComponent;
    class CameraComponent;
    class CameraControlComponent;
    class Input;
}

namespace Lizeral{
    class CameraControlSystem{

    public:
        CameraControlSystem();
        virtual ~CameraControlSystem();

        void Tick(float deltaTime, Registry& registry);

    private:
        void UpdateCameraDir(Input& input, TransformComponent& trans, CameraControlComponent& cameraController);

        void UpdateCameraPosForFree(Input& input, TransformComponent& trans, CameraControlComponent& cameraController, float deltaTime);
    };
}