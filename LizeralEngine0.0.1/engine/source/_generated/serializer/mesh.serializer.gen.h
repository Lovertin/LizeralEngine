#pragma once
#include "runtime\function\res_type\mesh\mesh.h"
#include "_generated\serializer\resource.serializer.gen.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Mesh& instance);
    template<>
    Mesh& PSerializer::read(const PJson& json_context, Mesh& instance);
}//namespace Lizeral

