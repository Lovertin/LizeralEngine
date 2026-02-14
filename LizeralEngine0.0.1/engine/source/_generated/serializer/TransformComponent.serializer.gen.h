#pragma once
#include "runtime\function\ecs\components\Transform\TransformComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const TransformComponent& instance);
    template<>
    TransformComponent& PSerializer::read(const PJson& json_context, TransformComponent& instance);
}//namespace Lizeral

