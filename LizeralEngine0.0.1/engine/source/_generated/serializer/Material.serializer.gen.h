#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\res_type\Material\Material.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Material& instance);
    template<>
    Material& PSerializer::read(const PJson& json_context, Material& instance);
}//namespace Lizeral

