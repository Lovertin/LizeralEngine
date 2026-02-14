#pragma once
#include "runtime\core\math\axis_aligned.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const AxisAlignedBox& instance);
    template<>
    AxisAlignedBox& PSerializer::read(const PJson& json_context, AxisAlignedBox& instance);
}//namespace Lizeral

