#pragma once
#include "runtime/function/ecs/components/component.h" 

// 定义脏标记位
const uint32_t PHYS_DIRTY_MASS = 1 << 1;
const uint32_t PHYS_DIRTY_FRICTION = 1 << 2;
const uint32_t PHYS_DIRTY_KINEMATIC = 1 << 3;

namespace Lizeral{
    REFLECTION_TYPE(RigidBodyComponent)
    CLASS(RigidBodyComponent:public Component, WhiteListFields){
        REFLECTION_BODY(RigidBodyComponent)
    public:
        // ========================================================
        // 必须在这里使用宏来生成代码！
        // ========================================================
        BEGIN_REFLECTION_UPDATED()
            REFLECTION_BIND_DIRTY("m_mass",          PHYS_DIRTY_MASS)
            REFLECTION_BIND_DIRTY("m_friction",      PHYS_DIRTY_FRICTION)
            REFLECTION_BIND_DIRTY("m_is_kinematic",  PHYS_DIRTY_KINEMATIC)
        END_REFLECTION_UPDATED()
        // ========================================================

        // Setter (运行时代码调用)
        void setMass(float mass) { 
            if (m_mass != mass) { m_mass = mass; setDirty(PHYS_DIRTY_MASS); } 
        }

        void setFriction(float friction){ 
            if(m_friction!=friction){m_friction=friction; setDirty(PHYS_DIRTY_FRICTION); }
        }

        void setKinematic(bool k) {
            if (m_is_kinematic != k) { m_is_kinematic = k; setDirty(PHYS_DIRTY_KINEMATIC); }
        }

        // Getter
        float getMass() const { return m_mass; }
        float getFriction() const { return m_friction; }
        bool isKinematic() const { return m_is_kinematic; }
        
        // System 接口
        void setRuntimeBody(void* body) { m_runtime_body = body; }
        void* getRuntimeBody() const { return m_runtime_body; }
    
    private:
        META(Enable)
        float m_mass { 1.0f };
        
        META(Enable)
        float m_friction { 0.5f };
        
        META(Enable)
        bool m_is_kinematic { false };

        void* m_runtime_body { nullptr }; 
    };
}