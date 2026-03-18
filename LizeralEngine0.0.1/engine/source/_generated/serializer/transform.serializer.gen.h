#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\transform.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Transform& instance);
    template<>
    inline Transform& PSerializer::read(const PJson& json_context, Transform& instance);
}//namespace Lizeral

