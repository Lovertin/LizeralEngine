#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Transform\TransformComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const TransformComponent& instance);
    template<>
    inline TransformComponent& PSerializer::read(const PJson& json_context, TransformComponent& instance);
}//namespace Lizeral

