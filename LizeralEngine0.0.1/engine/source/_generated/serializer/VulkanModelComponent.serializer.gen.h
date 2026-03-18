#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Model\VulkanModelComponent.h"
#include "_generated\serializer\component.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const VulkanModelComponent& instance);
    template<>
    VulkanModelComponent& PSerializer::read(const PJson& json_context, VulkanModelComponent& instance);
}//namespace Lizeral

