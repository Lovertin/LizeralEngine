#pragma once
#include "runtime\core\math\vector4.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Vector4& instance);
    template<>
    Vector4& PSerializer::read(const PJson& json_context, Vector4& instance);
}//namespace Lizeral

