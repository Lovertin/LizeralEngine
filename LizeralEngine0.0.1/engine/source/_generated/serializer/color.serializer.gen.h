#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\color\color.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Color& instance);
    template<>
    inline Color& PSerializer::read(const PJson& json_context, Color& instance);
}//namespace Lizeral

