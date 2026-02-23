#pragma once

// Bullet Headers
#include "btBulletDynamicsCommon.h"

// ECS & Components Headers
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/component.h" // 包含脏标记定义
#include "runtime/function/physics/phsicsEntityheaders.h"

// 前置声明，减少编译依赖
namespace Lizeral {
    class PhysicsScene;
    class RigidBodyComponent;
    class TransformComponent;
    class ColliderComponent;
    class PhysicsDebugDrawer;
}

namespace Lizeral {

    class PhysicsSystem {
    public:
        PhysicsSystem();
        virtual ~PhysicsSystem();

        // 初始化：绑定物理场景
        void Initialize(PhysicsScene* scene);

        // 核心循环：在 GameLoop 中每帧调用
        void Tick(float deltaTime, Registry& registry);

        // 清理
        void Shutdown();

        bool Raycast(const Ray& ray, RaycastHit& outHit, float maxDistance = 1000.0f);

        void DebugDrawWorld();

    private:
        PhysicsScene* m_scene { nullptr };

        PhysicsDebugDrawer* m_debugDrawer = nullptr;

    private:
        // --- 内部核心逻辑 ---

        // 1. 工厂方法：根据 ColliderComponent 创建 Bullet 形状
        btCollisionShape* CreateShape(const ColliderComponent& col);

        btCollisionShape* CreateShape(Entity entity, Registry& registry, const ColliderComponent& col);

        // 2. 创建刚体：将 ECS 组件转换为 Bullet 刚体并加入世界
        void CreateBody(Registry& registry,Entity entity, RigidBodyComponent& rb, const TransformComponent& t, const ColliderComponent& c);

        // 3. 数据同步 (Game -> Physics)：处理脏标记
        void SyncDirtyData(RigidBodyComponent& rb, const TransformComponent& t, ColliderComponent& c);

        // 4. 结果回写 (Physics -> Game)：将模拟结果写回组件
        void SyncTransformBack(const RigidBodyComponent& rb, TransformComponent& t,const ColliderComponent& c);

        // --- 数学辅助转换 ---
        btVector3 ToBtVector3(const class Vector3& v);
        btQuaternion ToBtQuaternion(const class Quaternion& q);
    };

} // namespace Lizeral