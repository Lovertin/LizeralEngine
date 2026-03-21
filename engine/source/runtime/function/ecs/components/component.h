#pragma once
#include <cstdint>
#include <cstring>
#include "runtime/core/meta/reflection/reflection.h"

namespace Lizeral{

    #define BEGIN_REFLECTION_UPDATED() \
    virtual void onReflectionUpdated(const char* field_name) override { \
        if (false) {} 

    #define REFLECTION_BIND_DIRTY(field_str, dirty_flag) \
        else if (std::strcmp(field_name, field_str) == 0) { \
            setDirty(dirty_flag); \
        }

    #define END_REFLECTION_UPDATED() \
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

        virtual void onReflectionUpdated(const char* field_name) {
            setDirty(DIRTY_GENERAL); 
        }

        bool isDirty() const {return m_dirty_flags!=0;}
        bool isDirty(uint32_t flag ) const { return (m_dirty_flags & flag) != 0; }
        void setDirty(uint32_t flag) { m_dirty_flags |= flag; }
        void clearDirty(uint32_t flag) { m_dirty_flags &= ~flag; }
    };
}