#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\vector2.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Vector2& instance);
    template<>
    inline Vector2& PSerializer::read(const PJson& json_context, Vector2& instance);
}//namespace Lizeral

