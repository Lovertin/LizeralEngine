#include "CameraControlSystem.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/input/input.h"
// #include <GLFW/glfw3.h> // 必须包含 GLFW 以使用 glfwSetInputMode

#include "runtime/core/math/vector2.h"

namespace Lizeral{

    CameraControlSystem::CameraControlSystem(){}

    CameraControlSystem::~CameraControlSystem(){}

    void CameraControlSystem::Tick(float deltaTime, Registry& registry) {
        auto& input = Input::GetInstance();

        // 1. 处理状态切换 (按下瞬间)
        if (input.GetMouseButtonDown(MouseButton::Right)) {
            // 【删除】：隐藏光标的代码已经移交给了 Qt 的 EngineViewportWidget 去做了！
            
            // 重置鼠标第一帧标记，防止出现巨大的镜头跳跃
            input.ResetMouse(); 
        }

        // 2. 处理持续的逻辑 (按住的过程中)
        // 只有按住右键，才允许更新相机的视角和位置
        if (input.GetMouseButton(MouseButton::Right)) {
            auto view = registry.view<TransformComponent, CameraControlComponent>();

            for (auto entity : view) {
                auto& trans = view.get<TransformComponent>(entity);
                auto& cameraControl = view.get<CameraControlComponent>(entity);

                UpdateCameraDir(input, trans, cameraControl);
                // 【修改】：把真实的 deltaTime 传进去
                UpdateCameraPosForFree(input, trans, cameraControl, deltaTime); 
            }
        }
    }

    void CameraControlSystem::UpdateCameraDir(Input& input,TransformComponent& trans,CameraControlComponent& cameraController){
        Vector2 delta = input.GetMouseDelta();

        if (delta.x == 0 && delta.y == 0) return;

        float newYaw = cameraController.getYaw() - delta.x * cameraController.getSensitivityX(); // 这里 m_sensitivityX 如果也是私有，也需要 getter
        float newPitch = cameraController.getPitch() - delta.y * cameraController.getSensitivityY();

        cameraController.setYaw(newYaw);
        cameraController.setPitch(newPitch);

        Quaternion qPitch = Quaternion::getQuaternionFromAngleAxis(Radian(Math::degreesToRadians(cameraController.getPitch())),Vector3(1.0f,0.0f,0.0f));
        Quaternion qYaw   = Quaternion::getQuaternionFromAngleAxis(Radian(Math::degreesToRadians(cameraController.getYaw())),Vector3(0.0f,1.0f,0.0f));

        trans.setRotation(qYaw * qPitch);
    }

    void CameraControlSystem::UpdateCameraPosForFree(Input& input, TransformComponent& trans, CameraControlComponent& cameraController, float dt){
        // 【修改】：使用传进来的真实 dt，而不是硬编码 0.016f
        float velocity = cameraController.getMoveSpeed() * dt;

        if (input.GetKey(Key::LEFT_SHIFT)) velocity = cameraController.getSpeedMutipier() * dt;

        Quaternion q = trans.getRotation();
        Vector3 forward = q.getForwardVector();
        Vector3 right   = q.getRightVector();

        Vector3 worldUp(0, 1, 0);
        Vector3 moveDir(0, 0, 0);

        // WASD - 沿着相机的局部坐标系移动
        if (input.GetKey(Key::W)) moveDir += forward; 
        if (input.GetKey(Key::S)) moveDir -= forward;
        if (input.GetKey(Key::D)) moveDir += right;
        if (input.GetKey(Key::A)) moveDir -= right;

        // QE - 沿着世界坐标系垂直升降
        if (input.GetKey(Key::Q)) moveDir += worldUp;
        if (input.GetKey(Key::E)) moveDir -= worldUp;

        if (moveDir.squaredLength() > 0) {
            moveDir.normalise();
            trans.setPosition(trans.getPosition() + moveDir * velocity);
        }
    }
}