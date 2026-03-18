#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Camera\CameraComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const CameraComponent& instance);
    template<>
    inline CameraComponent& PSerializer::read(const PJson& json_context, CameraComponent& instance);
}//namespace Lizeral

