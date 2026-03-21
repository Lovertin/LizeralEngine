#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\vector4.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Vector4& instance);
    template<>
    inline Vector4& PSerializer::read(const PJson& json_context, Vector4& instance);
}//namespace Lizeral

