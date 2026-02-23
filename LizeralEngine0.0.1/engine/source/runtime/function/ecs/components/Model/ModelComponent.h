#pragma once
#include <memory>
#include "runtime/function/ecs/components/component.h"
#include "runtime/function/res_type/Model/Model.h"

namespace Lizeral {
    REFLECTION_TYPE(ModelComponent)
    CLASS(ModelComponent:public Component,WhiteListFields) {
        REFLECTION_BODY(ModelComponent)
    public:
        std::shared_ptr<Model> m_Model;

        ModelComponent() = default;
        ModelComponent(std::shared_ptr<Model> model) : m_Model(model) {}
    };
}