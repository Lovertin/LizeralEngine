#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\test\parser_test.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const Test& instance);
    template<>
    inline Test& PSerializer::read(const PJson& json_context, Test& instance);
}//namespace Lizeral

