#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\component.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Component& instance);
    template<>
    inline Component& PSerializer::read(const PJson& json_context, Component& instance);
}//namespace Lizeral

