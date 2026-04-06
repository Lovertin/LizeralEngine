#pragma once
#include "hybrid_component_traits.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/RigidBody/RigidBodyComponent.h"
#include "runtime/function/ecs/components/Collider/ColliderComponent.h"

namespace Lizeral {

    template<>
    struct HybridStorageTrait<TransformComponent> {
        static constexpr HybridStorageKind storage = HybridStorageKind::Archetype;
    };

    template<>
    struct HybridStorageTrait<RigidBodyComponent> {
        static constexpr HybridStorageKind storage = HybridStorageKind::Archetype;
    };

    template<>
    struct HybridStorageTrait<ColliderComponent> {
        static constexpr HybridStorageKind storage = HybridStorageKind::Archetype;
    };

} // namespace Lizeral
