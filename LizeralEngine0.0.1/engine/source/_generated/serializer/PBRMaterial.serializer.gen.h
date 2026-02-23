#pragma once
#include "runtime\function\res_type\Material\PBRMaterial.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const PBRMaterial& instance);
    template<>
    PBRMaterial& PSerializer::read(const PJson& json_context, PBRMaterial& instance);
}//namespace Lizeral

