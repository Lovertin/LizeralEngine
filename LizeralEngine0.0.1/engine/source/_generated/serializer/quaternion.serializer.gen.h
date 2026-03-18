#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\quaternion.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Quaternion& instance);
    template<>
    Quaternion& PSerializer::read(const PJson& json_context, Quaternion& instance);
}//namespace Lizeral

