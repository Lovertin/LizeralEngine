#pragma once
#include "runtime\function\ecs\components\Model\ModelComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const ModelComponent& instance);
    template<>
    ModelComponent& PSerializer::read(const PJson& json_context, ModelComponent& instance);
}//namespace Lizeral

