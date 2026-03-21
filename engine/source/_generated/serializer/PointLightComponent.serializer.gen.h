#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Light\PointLightComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const PointLightComponent& instance);
    template<>
    inline PointLightComponent& PSerializer::read(const PJson& json_context, PointLightComponent& instance);
}//namespace Lizeral

