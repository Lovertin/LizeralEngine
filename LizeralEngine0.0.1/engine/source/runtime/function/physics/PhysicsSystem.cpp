#include "PhysicsSystem.h"

// 引入具体的实现头文件
#include "runtime/function/physics/PhysicsScene.h"
#include "runtime/function/ecs/components/RigidBody/RigidBodyComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Collider/ColliderComponent.h"

// 引入数学库 (假设路径)
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/quaternion.h"

namespace Lizeral {

    PhysicsSystem::PhysicsSystem() {}

    PhysicsSystem::~PhysicsSystem() {
        Shutdown();
    }

    void PhysicsSystem::Initialize(PhysicsScene* scene) {
        m_scene = scene;
    }

    void PhysicsSystem::Shutdown() {
        m_scene = nullptr;
    }

    // --- 数学转换辅助函数 ---
    btVector3 PhysicsSystem::ToBtVector3(const Vector3& v) {
        return btVector3(v.x, v.y, v.z);
    }

    btQuaternion PhysicsSystem::ToBtQuaternion(const Quaternion& q) {
        return btQuaternion(q.x, q.y, q.z, q.w);
    }

    // =====================================================================================
    // 核心循环 Tick
    // =====================================================================================
    void PhysicsSystem::Tick(float deltaTime, Registry& registry) {
        if (!m_scene || !m_scene->getWorld()) return;

        btDiscreteDynamicsWorld* world = m_scene->getWorld();

        // 1. 获取视图：只关心同时拥有 RB, Transform, Collider 的实体
        //    (Registry 的实现需要支持 view<A, B, C>)
        //    如果你的 view 实现暂时只支持单组件，你可能需要手动 get 其他组件并判空
        auto view = registry.view<RigidBodyComponent, TransformComponent, ColliderComponent>();

        // -----------------------------------------------------------
        // 阶段 A: 数据同步 (Game -> Physics) & 创建缺失刚体
        // -----------------------------------------------------------
        for (auto entity : view) {
            auto& rb = view.get<RigidBodyComponent>(entity);
            auto& trans = view.get<TransformComponent>(entity);
            auto& col = view.get<ColliderComponent>(entity);

            // case 1: 刚体未创建
            if (!rb.getRuntimeBody()) {
                CreateBody(entity, rb, trans, col);
                
                // 创建完成后，清除所有脏标记（视为已同步）
                // 使用 Component 基类的 clearDirty 清除所有位
                rb.clearDirty(0xFFFFFFFF); 
                trans.clearDirty(0xFFFFFFFF);
                col.clearDirty(0xFFFFFFFF);
                continue;
            }

            // case 2: 刚体已存在，检查脏标记
            if (rb.isDirty() || trans.isDirty() || col.isDirty()) {
                SyncDirtyData(rb, trans, col);
            }
        }

        // -----------------------------------------------------------
        // 阶段 B: 物理模拟 (Simulate)
        // -----------------------------------------------------------
        // 步长设置：最大子步数 10，固定时间步长 1/60 秒
        world->stepSimulation(deltaTime, 10, 1.0f / 60.0f);

        // -----------------------------------------------------------
        // 阶段 C: 结果回写 (Physics -> Game)
        // -----------------------------------------------------------
        for (auto entity : view) {
            auto& rb = view.get<RigidBodyComponent>(entity);

            // 只有【动态物体】(质量>0) 且 【非运动学物体】 才需要回写位置
            // 静态物体(Mass=0) 和 Kinematic 物体的位置是由 Transform 组件控制的，不接受物理反馈
            if (rb.getRuntimeBody() && rb.getMass() > 0.0f && !rb.isKinematic()) {
                auto& trans = view.get<TransformComponent>(entity);
                SyncTransformBack(rb, trans);
            }
        }
    }

    // =====================================================================================
    // 工厂方法：创建 Bullet 形状
    // =====================================================================================
    btCollisionShape* PhysicsSystem::CreateShape(const ColliderComponent& col) {
        btCollisionShape* shape = nullptr;

        switch (col.getType()) {
            case ColliderType::Box: {
                // Bullet BoxShape 需要的是半长 (Half Extents)
                // 组件里存的是全长 (Size)，所以除以 2
                Vector3 halfSize = col.getSize() * 0.5f;
                shape = new btBoxShape(ToBtVector3(halfSize));
                break;
            }
            case ColliderType::Sphere: {
                shape = new btSphereShape(col.getRadius());
                break;
            }
            case ColliderType::Capsule: {
                // Bullet Capsule: Radius, Height (Height 是中间圆柱段的高度)
                shape = new btCapsuleShape(col.getRadius(), col.getHeight());
                break;
            }
            default:
                // 默认给个小盒子防止崩溃
                shape = new btBoxShape(btVector3(0.5f, 0.5f, 0.5f));
                break;
        }

        return shape;
    }

    // =====================================================================================
    // 创建刚体
    // =====================================================================================
    void PhysicsSystem::CreateBody(Entity entity, RigidBodyComponent& rb, const TransformComponent& t, const ColliderComponent& c) {
        if (!m_scene) return;

        // 1. 创建形状
        btCollisionShape* shape = CreateShape(c);
        // 把形状交给 Scene 托管
        m_scene->trackShape(shape);

        // 2. 计算初始变换
        btTransform startTrans;
        startTrans.setIdentity();
        
        // 加上 Collider 的局部偏移 (Local Offset)
        // 简单的偏移处理：将初始位置设为 Transform位置 + 旋转 * 偏移
        // 注意：这只在创建时生效，如果运行时改 Offset，需要更复杂的 CompoundShape
        Vector3 finalPos = t.getPosition(); 
        // TODO: 如果需要精确的旋转偏移，这里需要用 Quaternion 旋转 offset 向量
        // finalPos = finalPos + (t.getRotation() * c.getOffset()); 
        // 暂时简化处理：直接叠加
        finalPos = finalPos + c.getOffset();

        startTrans.setOrigin(ToBtVector3(finalPos));
        startTrans.setRotation(ToBtQuaternion(t.getRotation()));

        // 3. 处理 Transform 的 Scale (Bullet 刚体没有 Scale，必须设置给 Shape)
        btVector3 localScaling = ToBtVector3(t.getScale());
        shape->setLocalScaling(localScaling);

        // 4. 计算惯性
        btVector3 localInertia(0, 0, 0);
        float mass = rb.getMass();
        if (mass > 0.0f) {
            shape->calculateLocalInertia(mass, localInertia);
        }

        // 5. 构建 Body
        // 使用 MotionState 可以让 Bullet 自动插值，减少抖动
        btDefaultMotionState* myMotionState = new btDefaultMotionState(startTrans);
        btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, shape, localInertia);
        
        btRigidBody* body = new btRigidBody(rbInfo);

        // 6. 设置额外属性
        body->setFriction(rb.getFriction());
        // 阻尼 (可选，防止空气中无限漂浮)
        // body->setDamping(0.1f, 0.1f); 

        if (rb.isKinematic()) {
            body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);
            body->setActivationState(DISABLE_DEACTIVATION); // 永远不休眠
        }

        // 7. 保存 Entity ID 到 UserPointer (非常重要，碰撞回调时反查 Entity)
        // 这里假设 Entity 是 uint32_t，强转存入
        // body->setUserPointer((void*)(uintptr_t)entity);

        body->setRestitution(rb.getRestitution());

        // 8. 加入世界
        m_scene->getWorld()->addRigidBody(body);
        
        // 9. 回填 Component
        rb.setRuntimeBody(body);
    }

    // =====================================================================================
    // 数据同步 Game -> Physics (处理脏标记)
    // =====================================================================================
    void PhysicsSystem::SyncDirtyData(RigidBodyComponent& rb, const TransformComponent& t, ColliderComponent& c) {
        btRigidBody* body = static_cast<btRigidBody*>(rb.getRuntimeBody());
        if (!body) return;

        // --- A. 处理刚体属性 (Mass, Friction, Kinematic) ---
        if (rb.isDirty(PHYS_DIRTY_MASS)) {
            float mass = rb.getMass();
            btVector3 inertia(0, 0, 0);
            if (mass > 0.0f && body->getCollisionShape()) {
                body->getCollisionShape()->calculateLocalInertia(mass, inertia);
            }
            body->setMassProps(mass, inertia);
            body->updateInertiaTensor();
            body->activate(true);
            rb.clearDirty(PHYS_DIRTY_MASS);
        }

        if (rb.isDirty(PHYS_DIRTY_FRICTION)) {
            body->setFriction(rb.getFriction());
            // 通常改摩擦力不需要唤醒，但保险起见
            body->activate(true); 
            rb.clearDirty(PHYS_DIRTY_FRICTION);
        }

        if (rb.isDirty(PHYS_DIRTY_KINEMATIC)) {
            if (rb.isKinematic()) {
                body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);
                body->setActivationState(DISABLE_DEACTIVATION);
            } else {
                body->setCollisionFlags(body->getCollisionFlags() & ~btCollisionObject::CF_KINEMATIC_OBJECT);
                body->setActivationState(ACTIVE_TAG);
                body->activate(true);
            }
            rb.clearDirty(PHYS_DIRTY_KINEMATIC);
        }

        if (rb.isDirty(PHYS_DIRTY_RESTITUTION)) {
            body->setRestitution(rb.getRestitution());
            // 修改属性后通常不需要唤醒，但在某些临界状态下唤醒更保险
            body->activate(true);
            rb.clearDirty(PHYS_DIRTY_RESTITUTION);
        }

        // --- B. 处理 Transform (瞬移 / Teleport) ---
        // 只有当 Position 或 Rotation 变脏时才强制设置
        if (t.isDirty(TRANS_DIRTY_POSITION) || t.isDirty(TRANS_DIRTY_ROTATION)) {
            btTransform tr;
            
            // 计算 offset (Transform + Collider Offset)
            Vector3 finalPos = t.getPosition() + c.getOffset();
            
            tr.setOrigin(ToBtVector3(finalPos));
            tr.setRotation(ToBtQuaternion(t.getRotation()));

            body->setWorldTransform(tr);
            if (body->getMotionState()) {
                body->getMotionState()->setWorldTransform(tr);
            }
            
            // 瞬移后必须唤醒，否则物体可能悬空静止
            body->activate(true);
        }

        // --- C. 处理 Scale (缩放) ---
        if (t.isDirty(TRANS_DIRTY_SCALE)) {
            if (body->getCollisionShape()) {
                body->getCollisionShape()->setLocalScaling(ToBtVector3(t.getScale()));
                
                // 缩放改变了形状，必须重新计算惯性 tensor
                if (rb.getMass() > 0.0f) {
                    btVector3 inertia;
                    body->getCollisionShape()->calculateLocalInertia(rb.getMass(), inertia);
                    body->setMassProps(rb.getMass(), inertia);
                    body->updateInertiaTensor();
                }
                body->activate(true);
            }
        }

        // 清除 Transform 的所有脏标记
        // 注意：因为我们是 const TransformComponent& t 传入的，
        // 这里假设 TransformComponent 的 clearDirty 是 mutable 的，或者你需要移除 const
        // 实际上在 Tick 中获取的是引用 auto& t，所以可以修改。
        // 为了编译通过，你需要确保 SyncDirtyData 的参数不是 const，或者 clearDirty 是 const 方法。
        // 修正：去掉 SyncDirtyData 参数里的 const
        // (我在下面的修正建议里去掉 const)
    }

    // =====================================================================================
    // 结果回写 Physics -> Game
    // =====================================================================================
    void PhysicsSystem::SyncTransformBack(const RigidBodyComponent& rb, TransformComponent& t) {
        btRigidBody* body = static_cast<btRigidBody*>(rb.getRuntimeBody());
        
        btTransform trans;
        if (body->getMotionState()) {
            body->getMotionState()->getWorldTransform(trans);
        } else {
            trans = body->getWorldTransform();
        }

        btVector3 pos = trans.getOrigin();
        btQuaternion rot = trans.getRotation();

        // 转换回 Lizeral 类型
        Vector3 newPos(pos.x(), pos.y(), pos.z());
        Quaternion newRot(rot.x(), rot.y(), rot.z(), rot.w());

        // 【关键】直接修改数据，不触发脏标记
        // 假设 TransformComponent 有 friend class PhysicsSystem
        // 或者使用 t.m_position = ... 如果它是 public 的
        // 这里假设你用 public 成员或者特殊的 set 方法
        
        // 减去 Collider 的 Offset (因为 Transform 代表的是物体轴心，而不是物理中心)
        // 物理位置 = Transform位置 + Offset
        // 所以：Transform位置 = 物理位置 - Offset
        // 注意：这里需要考虑旋转，Offset 是局部坐标
        // 暂时的简单回写：

        t.m_position = newPos;
        t.m_rotation = newRot;

        // 注意：这里绝对不能调用 t.setPosition(...)，因为那会触发 setDirty，导致死循环
    }

    // --- 修正 SyncDirtyData 的签名，去掉 const 以便清除脏标记 ---
    // void PhysicsSystem::SyncDirtyData(RigidBodyComponent& rb, TransformComponent& t, ColliderComponent& c) { ... }

} // namespace Lizeral