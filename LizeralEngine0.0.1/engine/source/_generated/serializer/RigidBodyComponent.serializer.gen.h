#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\RigidBody\RigidBodyComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const RigidBodyComponent& instance);
    template<>
    RigidBodyComponent& PSerializer::read(const PJson& json_context, RigidBodyComponent& instance);
}//namespace Lizeral

