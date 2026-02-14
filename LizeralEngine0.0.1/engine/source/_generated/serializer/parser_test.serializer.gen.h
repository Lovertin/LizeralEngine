#pragma once
#include "runtime\core\test\parser_test.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Test& instance);
    template<>
    Test& PSerializer::read(const PJson& json_context, Test& instance);
}//namespace Lizeral

