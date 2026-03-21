#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Light\DirectionalLightComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const DirectionLightComponent& instance);
    template<>
    inline DirectionLightComponent& PSerializer::read(const PJson& json_context, DirectionLightComponent& instance);
}//namespace Lizeral

