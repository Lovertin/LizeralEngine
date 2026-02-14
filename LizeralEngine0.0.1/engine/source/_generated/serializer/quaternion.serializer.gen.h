#pragma once
#include "runtime\core\math\quaternion.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Quaternion& instance);
    template<>
    Quaternion& PSerializer::read(const PJson& json_context, Quaternion& instance);
}//namespace Lizeral

