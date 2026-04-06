#pragma once

#include <btBulletDynamicsCommon.h>
#include <cstddef>
#include <unordered_map>
#include <vector>

#include "runtime/function/ecs/registry.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/quaternion.h"

namespace Lizeral {
    class PhysicsScene;
    class HybridRegistry;
    class RigidBodyComponent;
    class TransformComponent;
    class ColliderComponent;

    struct Ray;
    struct RaycastHit;

    struct DebugLineVertex {
        Vector3 position;
        Vector3 color;
    };

    struct PhysicsAccessProfile {
        std::size_t entity_count {0};
        double view_build_us {0.0};
        double iteration_us {0.0};
        double total_us {0.0};
    };

    struct PhysicsTickProfile {
        std::size_t entity_count {0};
        std::size_t created_body_count {0};
        std::size_t dirty_sync_count {0};
        std::size_t transform_write_back_count {0};
        double pre_step_us {0.0};
        double simulation_us {0.0};
        double write_back_us {0.0};
        double total_us {0.0};
    };

    struct PhysicsAccessComparison {
        PhysicsAccessProfile sparse_set;
        PhysicsAccessProfile hybrid_archetype;
    };

    enum class PhysicsHybridExecutionMode : uint8_t {
        View = 0,
        Chunked = 1
    };
}

namespace Lizeral {

    class PhysicsSystem {
    public:
        PhysicsSystem();
        virtual ~PhysicsSystem();

        void Initialize(PhysicsScene* scene);
        void Shutdown();
        void ResetPhysicsWorld();
        void ResetPhysicsWorld(Registry& registry);
        void ResetPhysicsWorld(HybridRegistry& registry);

        void SetHybridExecutionMode(PhysicsHybridExecutionMode mode) { m_hybridExecutionMode = mode; }
        PhysicsHybridExecutionMode GetHybridExecutionMode() const { return m_hybridExecutionMode; }

        void Tick(float deltaTime, Registry& registry);
        void Tick(float deltaTime, HybridRegistry& registry);
        void TickChunked(float deltaTime, HybridRegistry& registry);

        PhysicsTickProfile TickProfiled(float deltaTime, Registry& registry);
        PhysicsTickProfile TickProfiled(float deltaTime, HybridRegistry& registry);
        PhysicsTickProfile TickProfiledChunked(float deltaTime, HybridRegistry& registry);

        PhysicsAccessProfile ProfileComponentAccess(Registry& registry);
        PhysicsAccessProfile ProfileComponentAccess(HybridRegistry& registry);
        PhysicsAccessProfile ProfileChunkAccess(HybridRegistry& registry);
        PhysicsAccessComparison CompareComponentAccess(Registry& sparseRegistry, HybridRegistry& hybridRegistry);

        void DebugDrawWorld(Registry& registry);
        void DebugDrawWorld(HybridRegistry& registry);

        bool Raycast(const Ray& ray, RaycastHit& outHit, float maxDistance = 1000.0f);
        void Clear();

        PhysicsScene* getScene() { return m_scene; }

        void UpdateDebugLines(Registry& registry);
        void UpdateDebugLines(HybridRegistry& registry);

        const std::vector<DebugLineVertex>& GetDebugLines() const { return m_debugLines; }
        const PhysicsTickProfile& GetLastSparseTickProfile() const { return m_lastSparseTickProfile; }
        const PhysicsTickProfile& GetLastHybridTickProfile() const { return m_lastHybridTickProfile; }
        const PhysicsTickProfile& GetLastHybridChunkTickProfile() const { return m_lastHybridChunkTickProfile; }

    private:
        PhysicsScene* m_scene { nullptr };
        std::vector<DebugLineVertex> m_debugLines;
        PhysicsTickProfile m_lastSparseTickProfile;
        PhysicsTickProfile m_lastHybridTickProfile;
        PhysicsTickProfile m_lastHybridChunkTickProfile;
        std::unordered_map<Entity, btRigidBody*> m_trackedBodies;
        PhysicsHybridExecutionMode m_hybridExecutionMode { PhysicsHybridExecutionMode::Chunked };

    private:
        template<typename RegistryT>
        PhysicsTickProfile TickInternal(float deltaTime, RegistryT& registry, PhysicsTickProfile* profileSink);

        template<typename RegistryT>
        PhysicsAccessProfile ProfileComponentAccessInternal(RegistryT& registry);

        PhysicsTickProfile TickChunkedInternal(float deltaTime, HybridRegistry& registry, PhysicsTickProfile* profileSink);
        PhysicsAccessProfile ProfileChunkAccessInternal(HybridRegistry& registry);

        template<typename RegistryT>
        void UpdateDebugLinesInternal(RegistryT& registry);

        template<typename RegistryT>
        void DebugDrawWorldInternal(RegistryT& registry);

        template<typename RegistryT>
        void CleanupStaleBodiesInternal(RegistryT& registry);

        template<typename RegistryT>
        void ResetPhysicsWorldInternal(RegistryT& registry);

        void DestroyBody(btRigidBody* body);
        btCollisionShape* CreateShape(const ColliderComponent& col);
        void CreateBody(Entity entity, RigidBodyComponent& rb, const TransformComponent& t, const ColliderComponent& c);
        void SyncDirtyData(RigidBodyComponent& rb, TransformComponent& t, ColliderComponent& c);
        void SyncTransformBack(const RigidBodyComponent& rb, TransformComponent& t, const ColliderComponent& c);
        btVector3 ToBtVector3(const Vector3& v);
        btQuaternion ToBtQuaternion(const Quaternion& q);
    };

} // namespace Lizeral
