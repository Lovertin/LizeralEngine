#pragma once
#include <array>
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/components/component.h"

namespace Lizeral{
    REFLECTION_TYPE(NameComponent)
    CLASS(NameComponent:public Component,WhiteListFields){
        REFLECTION_BODY(NameComponent)

    public:
        NameComponent() = default;
        NameComponent(const std::string& name) : m_name(name) {}

        const std::string& getName() const { return m_name; }
        void setName(const std::string& name) { m_name = name; }
        
    private:
        META(Enable,UI:Headline)
        std::string m_name{"Empty Entity"};
    };
}