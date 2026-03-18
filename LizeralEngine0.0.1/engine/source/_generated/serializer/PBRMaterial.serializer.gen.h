#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\res_type\Material\PBRMaterial.h"
#include "_generated\serializer\Material.serializer.gen.h"

namespace Lizeral{
    template<>
    inline PJson PSerializer::write(const PBRMaterial& instance);
    template<>
    inline PBRMaterial& PSerializer::read(const PJson& json_context, PBRMaterial& instance);
}//namespace Lizeral

