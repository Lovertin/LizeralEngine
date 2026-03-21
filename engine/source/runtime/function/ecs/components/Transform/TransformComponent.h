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
            Vector3 forward = forwardDir.normalisedCopy();
            Vector3 up = upDir;

            if (std::abs(forward.dotProduct(up)) > 0.999f) {
                up = Vector3(0.0f, 0.0f, 1.0f);
            }

            Vector3 right = forward.crossProduct(up).normalisedCopy();
            Vector3 trueUp = right.crossProduct(forward).normalisedCopy();

            Quaternion newRot;
            newRot.fromAxes(right, trueUp, -forward);
            newRot.normalise();

            setRotation(newRot);
        }

    private:

        META(Enable,UI:Slider,Infinite)
        Vector3 m_position { Vector3::ZERO };

        META(Enable,UI:Slider,Range:-360~360)
        Quaternion m_rotation { Quaternion::IDENTITY };

        META(Enable)
        Vector3 m_scale { Vector3::UNIT_SCALE };

    };
}