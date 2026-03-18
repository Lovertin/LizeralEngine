#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\matrix4.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Matrix4x4_& instance);
    template<>
    inline Matrix4x4_& PSerializer::read(const PJson& json_context, Matrix4x4_& instance);
}//namespace Lizeral

