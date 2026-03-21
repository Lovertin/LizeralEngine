#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\res_type\texture\TextureCube.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const TextureCube& instance);
    template<>
    inline TextureCube& PSerializer::read(const PJson& json_context, TextureCube& instance);
}//namespace Lizeral

