#pragma once
#include "runtime\function\ecs\components\component.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Component& instance);
    template<>
    Component& PSerializer::read(const PJson& json_context, Component& instance);
}//namespace Lizeral

