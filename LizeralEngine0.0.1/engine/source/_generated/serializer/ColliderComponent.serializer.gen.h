#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Collider\ColliderComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const ColliderComponent& instance);
    template<>
    ColliderComponent& PSerializer::read(const PJson& json_context, ColliderComponent& instance);
}//namespace Lizeral

