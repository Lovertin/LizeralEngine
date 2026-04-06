#pragma once

namespace Lizeral {

    enum class HybridStorageKind {
        SparseSet,
        Archetype
    };

    template<typename T>
    struct HybridStorageTrait {
        static constexpr HybridStorageKind storage = HybridStorageKind::SparseSet;
    };

    template<typename T>
    inline constexpr bool is_hybrid_archetype_component_v =
        HybridStorageTrait<T>::storage == HybridStorageKind::Archetype;

} // namespace Lizeral
