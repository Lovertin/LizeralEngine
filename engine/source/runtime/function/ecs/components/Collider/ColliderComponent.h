#pragma once
#include "runtime/function/ecs/components/component.h"
#include "runtime/core/math/vector3.h"

const uint32_t COLLIDER_DIRTY_TYPE   = 1 << 1; 
const uint32_t COLLIDER_DIRTY_SIZE   = 1 << 2; 
const uint32_t COLLIDER_DIRTY_OFFSET = 1 << 3; 

namespace Lizeral {

    enum class ColliderType {
        Box,
        Sphere,
        Capsule,
        ConvexHull
    };

    REFLECTION_TYPE(ColliderComponent)
    CLASS(ColliderComponent : public Component, WhiteListFields) {
        REFLECTION_BODY(ColliderComponent)

    public:
        BEGIN_REFLECTION_UPDATED()
            REFLECTION_BIND_DIRTY("m_type",     COLLIDER_DIRTY_TYPE)
            REFLECTION_BIND_DIRTY("m_size",     COLLIDER_DIRTY_SIZE)
            REFLECTION_BIND_DIRTY("m_radius",   COLLIDER_DIRTY_SIZE)
            REFLECTION_BIND_DIRTY("m_height",   COLLIDER_DIRTY_SIZE)
            REFLECTION_BIND_DIRTY("m_offset",   COLLIDER_DIRTY_OFFSET)
        END_REFLECTION_UPDATED()

        // --- Setters ---
        void setType(ColliderType type) {
            if (m_type != type) { m_type = type; setDirty(COLLIDER_DIRTY_TYPE); }
        }

        void setSize(const Vector3& size) { // Box
            if (m_size != size) { m_size = size; setDirty(COLLIDER_DIRTY_SIZE); }
        }

        void setRadius(float radius) { // Sphere, Capsule
            if (m_radius != radius) { m_radius = radius; setDirty(COLLIDER_DIRTY_SIZE); }
        }

        void setHeight(float height) { // Capsule
            if (m_height != height) { m_height = height; setDirty(COLLIDER_DIRTY_SIZE); }
        }

        void setOffset(const Vector3& offset) {
            if (m_offset != offset) { m_offset = offset; setDirty(COLLIDER_DIRTY_OFFSET); }
        }

        void setDebug(const bool isDebug){ m_ShowDebug = isDebug; }

        // --- Getters ---
        ColliderType getType() const { return m_type; }
        const Vector3& getSize() const { return m_size; }
        float getRadius() const { return m_radius; }
        float getHeight() const { return m_height; }
        const Vector3& getOffset() const { return m_offset; }
        bool getShowDebug() const {return m_ShowDebug;}

    private:
        
        META(Enable)
        ColliderType m_type { ColliderType::Box };

        META(Enable)
        Vector3 m_size { 1.0f, 1.0f, 1.0f };

        META(Enable)
        float m_radius { 0.5f };

        META(Enable)
        float m_height { 1.0f };

        META(Enable)
        Vector3 m_offset { 0.0f, 0.0f, 0.0f };

        META(Enable)
        bool m_ShowDebug {false};
    };
}