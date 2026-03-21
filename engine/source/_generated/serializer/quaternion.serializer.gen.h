#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\quaternion.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Quaternion& instance);
    template<>
    inline Quaternion& PSerializer::read(const PJson& json_context, Quaternion& instance);
}//namespace Lizeral

