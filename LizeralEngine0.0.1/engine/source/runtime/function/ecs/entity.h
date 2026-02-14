#pragma once
#include <cstdint>

namespace Lizeral {

    // 简单起见，使用 uint32_t 作为 ID
    // 高级实现通常会用高位做 Version 检查以防止 ID 复用 bug，这里简化处理
    using Entity = std::uint32_t;

    // 定义一个无效的 Entity ID
    constexpr Entity null_entity = static_cast<Entity>(-1);

} // namespace Lizeral