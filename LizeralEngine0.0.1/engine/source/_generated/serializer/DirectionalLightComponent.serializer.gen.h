#pragma once
#include "runtime\function\ecs\components\Light\DirectionalLightComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const DirectionLightComponent& instance);
    template<>
    DirectionLightComponent& PSerializer::read(const PJson& json_context, DirectionLightComponent& instance);
}//namespace Lizeral

