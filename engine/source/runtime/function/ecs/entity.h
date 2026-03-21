#pragma once
#include <cstdint>

namespace Lizeral {

    using Entity = std::uint32_t;

    // define Entity ID
    constexpr Entity null_entity = static_cast<Entity>(-1);

} // namespace Lizeral