#pragma once
#include "runtime\resource\resource.h"

namespace Lizeral{
    template<>
    PJson PSerializer::write(const Resource& instance);
    template<>
    Resource& PSerializer::read(const PJson& json_context, Resource& instance);
}//namespace Lizeral

