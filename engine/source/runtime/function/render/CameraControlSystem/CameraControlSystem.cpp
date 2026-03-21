#include "CameraControlSystem.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/input/input.h"
// #include <GLFW/glfw3.h> 

#include "runtime/core/math/vector2.h"

namespace Lizeral{

    CameraControlSystem::CameraControlSystem(){}

    CameraControlSystem::~CameraControlSystem(){}

    void CameraControlSystem::Tick(float deltaTime, Registry& registry) {
        auto& input = Input::GetInstance();

        if (input.GetMouseButtonDown(MouseButton::Right)) {

            input.ResetMouse(); 
        }

        if (input.GetMouseButton(MouseButton::Right)) {
            auto view = registry.view<TransformComponent, CameraControlComponent>();

            for (auto entity : view) {
                auto& trans = view.get<TransformComponent>(entity);
                auto& cameraControl = view.get<CameraControlComponent>(entity);

                UpdateCameraDir(input, trans, cameraControl);

                UpdateCameraPosForFree(input, trans, cameraControl, deltaTime); 
            }
        }
    }

    void CameraControlSystem::UpdateCameraDir(Input& input,TransformComponent& trans,CameraControlComponent& cameraController){
        Vector2 delta = input.GetMouseDelta();

        if (delta.x == 0 && delta.y == 0) return;

        float newYaw = cameraController.getYaw() - delta.x * cameraController.getSensitivityX();
        float newPitch = cameraController.getPitch() - delta.y * cameraController.getSensitivityY();

        cameraController.setYaw(newYaw);
        cameraController.setPitch(newPitch);

        Quaternion qPitch = Quaternion::getQuaternionFromAngleAxis(Radian(Math::degreesToRadians(cameraController.getPitch())),Vector3(1.0f,0.0f,0.0f));
        Quaternion qYaw   = Quaternion::getQuaternionFromAngleAxis(Radian(Math::degreesToRadians(cameraController.getYaw())),Vector3(0.0f,1.0f,0.0f));

        trans.setRotation(qYaw * qPitch);
    }

    void CameraControlSystem::UpdateCameraPosForFree(Input& input, TransformComponent& trans, CameraControlComponent& cameraController, float dt){
        float velocity = cameraController.getMoveSpeed() * dt;

        if (input.GetKey(Key::LEFT_SHIFT)) velocity = cameraController.getSpeedMutipier() * dt;

        Quaternion q = trans.getRotation();
        Vector3 forward = q.getForwardVector();
        Vector3 right   = q.getRightVector();

        Vector3 worldUp(0, 1, 0);
        Vector3 moveDir(0, 0, 0);

        // WASD forward
        if (input.GetKey(Key::W)) moveDir += forward; 
        if (input.GetKey(Key::S)) moveDir -= forward;
        if (input.GetKey(Key::D)) moveDir += right;
        if (input.GetKey(Key::A)) moveDir -= right;

        // QE up & down
        if (input.GetKey(Key::Q)) moveDir += worldUp;
        if (input.GetKey(Key::E)) moveDir -= worldUp;

        if (moveDir.squaredLength() > 0) {
            moveDir.normalise();
            trans.setPosition(trans.getPosition() + moveDir * velocity);
        }
    }
}