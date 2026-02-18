#include "CameraControlSystem.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"

#include "runtime/core/math/vector2.h"

namespace Lizeral{

    CameraControlSystem::CameraControlSystem(){}

    CameraControlSystem::~CameraControlSystem(){}

    void CameraControlSystem::Tick(Registry& registry){
        auto view=registry.view<TransformComponent,CameraControlComponent>();

        for(auto entity:view){
            auto& input=Input::GetInstance();
            auto& trans=view.get<TransformComponent>(entity);
            auto& cameraControl=view.get<CameraControlComponent>(entity);

            //每一帧更新方向和位置
            UpdateCameraDir(input,trans,cameraControl);

            UpdateCameraPosForFree(input,trans,cameraControl);
        }
    }

    void CameraControlSystem::UpdateCameraDir(Input& input,TransformComponent& trans,CameraControlComponent& cameraController){
        Vector2 delta = input.GetMouseDelta();

        if (delta.x == 0 && delta.y == 0) return;

        float newYaw = cameraController.getYaw() - delta.x * cameraController.getSensitivityX(); // 这里 m_sensitivityX 如果也是私有，也需要 getter
        float newPitch = cameraController.getPitch() - delta.y * cameraController.getSensitivityY();

        cameraController.setYaw(newYaw);
        cameraController.setPitch(newPitch);

        Quaternion qPitch = Quaternion::getQuaternionFromAngleAxis(Radian(Math::degreesToRadians(cameraController.getPitch())),Vector3(-1.0f,0.0f,0.0f));
        Quaternion qYaw   = Quaternion::getQuaternionFromAngleAxis(Radian(Math::degreesToRadians(cameraController.getYaw())),Vector3(0.0f,1.0f,0.0f));

        trans.setRotation(qYaw * qPitch);
    }

    void CameraControlSystem::UpdateCameraPosForFree(Input& input, TransformComponent& trans, CameraControlComponent& cameraController){
        float dt = 0.016f; // 这里应该从 Tick 传入 deltaTime
        float velocity = cameraController.getMoveSpeed() * dt;

        if (input.GetKey(Key::LEFT_SHIFT)) velocity = cameraController.getSpeedMutipier() * dt;

        Quaternion q = trans.getRotation();
        Vector3 forward = q.getForwardVector();
        Vector3 right   = q.getRightVector();

        Vector3 worldUp(0, 1, 0);

        Vector3 moveDir(0, 0, 0);

        // WASD - 沿着相机的局部坐标系移动
        if (input.GetKey(Key::W)) moveDir -= forward; //虽然反直觉，但是opengl中z的负半轴才是正方向
        if (input.GetKey(Key::S)) moveDir += forward;
        if (input.GetKey(Key::D)) moveDir += right;
        if (input.GetKey(Key::A)) moveDir -= right;

        // QE - 沿着世界坐标系垂直升降 (UE 风格)
        if (input.GetKey(Key::Q)) moveDir += worldUp;
        if (input.GetKey(Key::E)) moveDir -= worldUp;

        if (moveDir.squaredLength() > 0) {
            moveDir.normalise();
            trans.setPosition(trans.getPosition() + moveDir * velocity);
        }
    }
}