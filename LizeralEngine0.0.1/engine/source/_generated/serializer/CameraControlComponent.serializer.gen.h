#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Camera\CameraControlComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const CameraControlComponent& instance);
    template<>
    inline CameraControlComponent& PSerializer::read(const PJson& json_context, CameraControlComponent& instance);
}//namespace Lizeral

