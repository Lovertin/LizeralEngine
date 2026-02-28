#pragma once
#include "runtime\function\res_type\Model\Mesh.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Mesh& instance);
    template<>
    Mesh& PSerializer::read(const PJson& json_context, Mesh& instance);
}//namespace Lizeral

