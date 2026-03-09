#pragma once
#include "runtime/function/ecs/components/component.h"
#include "runtime/core/math/vector3.h"

// --- 定义脏标记 (从 1 开始，不与 RB 混用) ---
const uint32_t COLLIDER_DIRTY_TYPE   = 1 << 1; // 类型变了 (需要重建 Shape)
const uint32_t COLLIDER_DIRTY_SIZE   = 1 << 2; // 尺寸变了 (Box长宽/球半径/胶囊高度)
const uint32_t COLLIDER_DIRTY_OFFSET = 1 << 3; // 偏移变了

namespace Lizeral {

    // 1. 定义支持的形状枚举
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

        void setSize(const Vector3& size) { // Box 用
            if (m_size != size) { m_size = size; setDirty(COLLIDER_DIRTY_SIZE); }
        }

        void setRadius(float radius) { // Sphere, Capsule 用
            if (m_radius != radius) { m_radius = radius; setDirty(COLLIDER_DIRTY_SIZE); }
        }

        void setHeight(float height) { // Capsule 用
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
        // 2. 核心数据区
        
        META(Enable)
        ColliderType m_type { ColliderType::Box };

        META(Enable)
        Vector3 m_size { 1.0f, 1.0f, 1.0f };

        // Sphere 和 Capsule 的半径
        META(Enable)
        float m_radius { 0.5f };

        // Capsule 的高度 (通常指中间圆柱段的高度)
        META(Enable)
        float m_height { 1.0f };

        // 局部中心偏移 (相对于 Transform 的位置)
        META(Enable)
        Vector3 m_offset { 0.0f, 0.0f, 0.0f };

        META(Enable)
        bool m_ShowDebug {false};
    };
}