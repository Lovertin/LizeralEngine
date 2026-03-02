#pragma once
#include "runtime\function\res_type\texture\TextureCube.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const TextureCube& instance);
    template<>
    TextureCube& PSerializer::read(const PJson& json_context, TextureCube& instance);
}//namespace Lizeral

