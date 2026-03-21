#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\axis_aligned.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const AxisAlignedBox& instance);
    template<>
    inline AxisAlignedBox& PSerializer::read(const PJson& json_context, AxisAlignedBox& instance);
}//namespace Lizeral

