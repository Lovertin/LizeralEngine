#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Name\NameComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const NameComponent& instance);
    template<>
    inline NameComponent& PSerializer::read(const PJson& json_context, NameComponent& instance);
}//namespace Lizeral

