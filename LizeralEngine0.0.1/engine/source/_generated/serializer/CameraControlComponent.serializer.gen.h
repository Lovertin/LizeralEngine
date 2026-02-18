#pragma once
#include "runtime\function\ecs\components\Camera\CameraControlComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const CameraControlComponent& instance);
    template<>
    CameraControlComponent& PSerializer::read(const PJson& json_context, CameraControlComponent& instance);
}//namespace Lizeral

