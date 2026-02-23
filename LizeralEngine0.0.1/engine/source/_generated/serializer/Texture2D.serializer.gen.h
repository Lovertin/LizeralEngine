#pragma once
#include "runtime\function\res_type\texture\Texture2D.h"
#include "_generated\serializer\resource.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Texture2D& instance);
    template<>
    Texture2D& PSerializer::read(const PJson& json_context, Texture2D& instance);
}//namespace Lizeral

