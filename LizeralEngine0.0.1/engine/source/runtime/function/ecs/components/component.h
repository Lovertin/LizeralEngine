#pragma once
#include <cstdint>
#include <cstring>
#include "runtime/core/meta/reflection/reflection.h"

namespace Lizeral{

    // --- 宏定义放在这里没问题 ---
    #define BEGIN_REFLECTION_UPDATED() \
    virtual void onReflectionUpdated(const char* field_name) override { \
        if (false) {} /* 占位符，为了让后面的 else if 语法合法 */

    #define REFLECTION_BIND_DIRTY(field_str, dirty_flag) \
        else if (std::strcmp(field_name, field_str) == 0) { \
            setDirty(dirty_flag); \
        }

    #define END_REFLECTION_UPDATED() \
        /* 调用基类处理剩余逻辑 (链式调用) */ \
        Component::onReflectionUpdated(field_name); \
    }

    REFLECTION_TYPE(Component)
    CLASS(Component, WhiteListFields){
        REFLECTION_BODY(Component)
    protected:
        uint32_t m_dirty_flags {0};

    public:
        static const uint32_t DIRTY_NONE = 0;
        static const uint32_t DIRTY_GENERAL = 1 << 0;
        
        //  基类必须定义这个虚函数，否则子类无法 override
        virtual void onReflectionUpdated(const char* field_name) {
            // 基类默认把 DIRTY_GENERAL 设为 true，作为保底
            setDirty(DIRTY_GENERAL); 
        }

        bool isDirty(uint32_t flag = DIRTY_GENERAL) const { return (m_dirty_flags & flag) != 0; }
        void setDirty(uint32_t flag) { m_dirty_flags |= flag; }
        void clearDirty(uint32_t flag) { m_dirty_flags &= ~flag; }
    };
}