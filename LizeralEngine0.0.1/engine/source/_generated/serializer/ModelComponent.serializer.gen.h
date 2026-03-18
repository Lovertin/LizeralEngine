#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Model\ModelComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const ModelComponent& instance);
    template<>
    inline ModelComponent& PSerializer::read(const PJson& json_context, ModelComponent& instance);
}//namespace Lizeral

