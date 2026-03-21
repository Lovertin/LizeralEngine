#pragma once

// Bullet Headers
#include <btBulletDynamicsCommon.h>

// ECS & Math Headers
#include "runtime/function/ecs/registry.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/quaternion.h"

// 前置声明，极大加快编译速度，切断循环依赖
namespace Lizeral {
    class PhysicsScene;
    class RigidBodyComponent;
    class TransformComponent;
    class ColliderComponent;
    
    struct Ray;
    struct RaycastHit;

    struct DebugLineVertex {
        Vector3 position;
        Vector3 color;
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

        void Tick(float deltaTime, Registry& registry);
        void DebugDrawWorld(Registry& registry);
        bool Raycast(const Ray& ray, RaycastHit& outHit, float maxDistance = 1000.0f);
        void Clear();

        PhysicsScene* getScene() { return m_scene; }

        void UpdateDebugLines(Registry& registry);

        const std::vector<DebugLineVertex>& GetDebugLines() const { return m_debugLines; }


    private:
        PhysicsScene* m_scene { nullptr };

        std::vector<DebugLineVertex> m_debugLines;

    private:

        btCollisionShape* CreateShape(const ColliderComponent& col);
        btCollisionShape* CreateShape(Entity entity, Registry& registry, const ColliderComponent& col);

        void CreateBody(Registry& registry, Entity entity, RigidBodyComponent& rb, const TransformComponent& t, const ColliderComponent& c);

        void SyncDirtyData(RigidBodyComponent& rb, TransformComponent& t, ColliderComponent& c);

        void SyncTransformBack(const RigidBodyComponent& rb, TransformComponent& t, const ColliderComponent& c);

        btVector3 ToBtVector3(const Vector3& v);
        btQuaternion ToBtQuaternion(const Quaternion& q);
    };

} // namespace Lizeral