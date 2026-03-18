#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\color\color.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Color& instance);
    template<>
    Color& PSerializer::read(const PJson& json_context, Color& instance);
}//namespace Lizeral

