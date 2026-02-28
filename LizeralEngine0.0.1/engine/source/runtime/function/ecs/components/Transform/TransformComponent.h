#pragma once
#include "runtime/function/ecs/components/component.h" 
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/quaternion.h"
#include "runtime/core/math/matrix4.h"

namespace Lizeral{
    const uint32_t TRANS_DIRTY_POSITION = 1 << 1;
    const uint32_t TRANS_DIRTY_ROTATION = 1 << 2;
    const uint32_t TRANS_DIRTY_SCALE    = 1 << 3;

    REFLECTION_TYPE(TransformComponent)
    CLASS(TransformComponent:public Component,WhiteListFields){
        REFLECTION_BODY(TransformComponent)
    public:
        BEGIN_REFLECTION_UPDATED()
            REFLECTION_BIND_DIRTY("m_position", TRANS_DIRTY_POSITION)
            REFLECTION_BIND_DIRTY("m_rotation", TRANS_DIRTY_ROTATION)
            REFLECTION_BIND_DIRTY("m_scale",    TRANS_DIRTY_SCALE)
        END_REFLECTION_UPDATED()

        friend class PhysicsSystem;

        // Setters
        void setPosition(const Vector3& pos) {
            if (m_position != pos) { m_position = pos; setDirty(TRANS_DIRTY_POSITION); }
        }
        
        void setRotation(const Quaternion& rot) {
            if (m_rotation != rot) { m_rotation = rot; setDirty(TRANS_DIRTY_ROTATION); }
        }

        void setScale(const Vector3& scale) {
            if (m_scale != scale) { m_scale = scale; setDirty(TRANS_DIRTY_SCALE); }
        }

        // Getters ...
        const Vector3& getPosition() const { return m_position; }
        const Quaternion& getRotation() const { return m_rotation; }
        const Vector3& getScale() const { return m_scale; }

        // 获取变换矩阵（用于渲染）
        Matrix4x4 getMatrix() const {
             Matrix4x4 mat;
             mat.makeTransform(m_position, m_scale, m_rotation);
             return mat;
        }

        Vector3 getForward() const{
            float x = 2.0f * (m_rotation.x * m_rotation.z + m_rotation.w * m_rotation.y);
            float y = 2.0f * (m_rotation.y * m_rotation.z - m_rotation.w * m_rotation.x);
            float z = 1.0f - 2.0f * (m_rotation.x * m_rotation.x + m_rotation.y * m_rotation.y);
            return Vector3(-x, -y, -z).normalisedCopy();
        }

        Vector3 getUp() const {
            return m_rotation * Vector3(0.0f, 1.0f, 0.0f);
        }

        Vector3 getRight() const {
            return m_rotation * Vector3(1.0f, 0.0f, 0.0f);
        }

        void setForward(const Vector3& forwardDir, const Vector3& upDir = Vector3(0.0f, 1.0f, 0.0f)) {
            // 1. 获取标准化的前向向量
            Vector3 forward = forwardDir.normalisedCopy();
            Vector3 up = upDir;

            // 2. 防止万向节死锁：如果光线是垂直往下（或往上）打的，
            // 此时 forward 和 up 几乎平行，会导致叉乘结果为 0。我们需要临时换一个 up 向量。
            // (假设你的 Vector3 有 dotProduct 方法，如果没有可以写成 x*x + y*y + z*z)
            if (std::abs(forward.dotProduct(up)) > 0.999f) {
                up = Vector3(0.0f, 0.0f, 1.0f);
            }

            // 3. 叉乘构建正交基底 (右向量和真·上向量)
            Vector3 right = forward.crossProduct(up).normalisedCopy();
            Vector3 trueUp = right.crossProduct(forward).normalisedCopy();

            // 4. 构建四元数
            // 根据你 getForward 的实现，你的局部 Forward 对应的是 -Z 轴
            Quaternion newRot;
            newRot.fromAxes(right, trueUp, -forward);
            newRot.normalise();

            // 5. 调用已有的 setRotation
            // 【极其关键】：这里复用了 setRotation，所以它会自动为你打上 TRANS_DIRTY_ROTATION 的脏标记！
            setRotation(newRot);
        }

    private:
        META(Enable)
        Vector3 m_position { Vector3::ZERO };

        META(Enable)
        Quaternion m_rotation { Quaternion::IDENTITY };

        META(Enable)
        Vector3 m_scale { Vector3::UNIT_SCALE };

    };
}