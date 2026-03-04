#pragma once
#include "runtime\function\ecs\components\Name\NameComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const NameComponent& instance);
    template<>
    NameComponent& PSerializer::read(const PJson& json_context, NameComponent& instance);
}//namespace Lizeral

