#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\vector3.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Vector3& instance);
    template<>
    inline Vector3& PSerializer::read(const PJson& json_context, Vector3& instance);
}//namespace Lizeral

