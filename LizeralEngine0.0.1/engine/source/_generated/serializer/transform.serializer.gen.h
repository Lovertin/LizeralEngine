#pragma once
#include "runtime\core\math\transform.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Transform& instance);
    template<>
    Transform& PSerializer::read(const PJson& json_context, Transform& instance);
}//namespace Lizeral

