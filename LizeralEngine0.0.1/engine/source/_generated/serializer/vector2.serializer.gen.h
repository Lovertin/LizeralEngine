#pragma once
#include "runtime\core\math\vector2.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Vector2& instance);
    template<>
    Vector2& PSerializer::read(const PJson& json_context, Vector2& instance);
}//namespace Lizeral

