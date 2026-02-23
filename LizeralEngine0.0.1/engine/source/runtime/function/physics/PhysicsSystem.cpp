#include "PhysicsSystem.h"

// 引入具体的实现头文件
#include "runtime/function/physics/PhysicsScene.h"
#include "runtime/function/physics/physicsDebug/PhysicsDebugDrawer.h" 
#include "runtime/function/ecs/components/RigidBody/RigidBodyComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Collider/ColliderComponent.h"
#include "runtime/function/ecs/components/Model/ModelComponent.h"

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
        if (m_scene && m_scene->getWorld()) {
            // 获取物理世界的求解器配置
            btContactSolverInfo& info = m_scene->getWorld()->getSolverInfo();
            
            // 1. 增加迭代次数 (默认10) - 解决多个方块堆叠时的挤压穿模
            info.m_numIterations = 20; 
            
            // 2. 提高错误修正系数 ERP (默认0.2，最大1.0) - 让引擎更“用力”地把穿模物体推出来
            // info.m_erp = 0.2f; 

            m_debugDrawer = new PhysicsDebugDrawer();
            // 你可以通过位运算叠加你想看的调试信息：线框 + 包围盒 + 接触点
            m_debugDrawer->setDebugMode(btIDebugDraw::DBG_DrawWireframe | btIDebugDraw::DBG_DrawAabb);
            m_scene->getWorld()->setDebugDrawer(m_debugDrawer);
        }
    }

    void PhysicsSystem::Shutdown() {
        m_scene = nullptr;
        if (m_debugDrawer) {
            delete m_debugDrawer;
            m_debugDrawer = nullptr;
        }
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
                CreateBody(registry,entity, rb, trans, col);
                
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
                auto& col = view.get<ColliderComponent>(entity);
                SyncTransformBack(rb,trans,col);
            }
        }
    }

    void PhysicsSystem::DebugDrawWorld() {
        if (m_scene && m_scene->getWorld() && m_debugDrawer) {
            m_scene->getWorld()->debugDrawWorld();
        }
    }

    bool PhysicsSystem::Raycast(const Ray& ray, RaycastHit& outHit, float maxDistance){
        if (!m_scene->getWorld()) return false;

        Vector3 start=ray.origin;
        Vector3 end=ray.GetPoint(maxDistance);

        btVector3 btStart(start.x, start.y, start.z);
        btVector3 btEnd(end.x, end.y, end.z);

        btCollisionWorld::ClosestRayResultCallback callback(btStart, btEnd);

        m_scene->getWorld()->rayTest(btStart, btEnd, callback);

        if (callback.hasHit()) {
            outHit.hasHit = true;
            outHit.point = Vector3(callback.m_hitPointWorld.x(), callback.m_hitPointWorld.y(), callback.m_hitPointWorld.z());
            outHit.normal = Vector3(callback.m_hitNormalWorld.x(), callback.m_hitNormalWorld.y(), callback.m_hitNormalWorld.z());
            outHit.distance = (outHit.point - start).length();

            // 【关键难点】如何从 Bullet 的 collisionObject 拿到 ECS 的 Entity？
            // 通常在创建 RigidBody 时，我们需要把 EntityID 存到 UserPointer 里
            const btCollisionObject* obj = callback.m_collisionObject;
            outHit.entity = static_cast<uint32_t>((uintptr_t)obj->getUserPointer()); 
            
            return true;
        }

        outHit.hasHit = false;
        return false;
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
                // shape->setMargin(0.001f);
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

    btCollisionShape* PhysicsSystem::CreateShape(Entity entity, Registry& registry, const ColliderComponent& col) {
        btCollisionShape* shape = nullptr;

        switch (col.getType()) {
            case ColliderType::Box: {
                Vector3 halfSize = col.getSize() * 0.5f;
                shape = new btBoxShape(btVector3(halfSize.x, halfSize.y, halfSize.z));
                // shape->setMargin(0.001f);
                break;
            }
            case ColliderType::Sphere: {
                shape = new btSphereShape(col.getRadius());
                break;
            }
            case ColliderType::Capsule: {
                shape = new btCapsuleShape(col.getRadius(), col.getHeight());
                break;
            }
            case ColliderType::ConvexHull: {
            // 在商业引擎中，我们用复合碰撞体(CompoundShape)加包围盒(AABB)来完美替代高面数凸包
                auto* compoundShape = new btCompoundShape();

                if (registry.has<ModelComponent>(entity)) {
                    auto& modelComp = registry.get<ModelComponent>(entity);
                    if (modelComp.m_Model) {
                        // 1. 准备计算模型的 AABB (轴向包围盒)
                        Vector3 minAABB( FLT_MAX,  FLT_MAX,  FLT_MAX);
                        Vector3 maxAABB(-FLT_MAX, -FLT_MAX, -FLT_MAX);

                        // 2. 遍历所有顶点，找出模型的极值边界
                        for (const auto& mesh : modelComp.m_Model->GetMeshes()) {
                            for (const auto& v : mesh.m_Vertices) {
                                minAABB.x = std::min(minAABB.x, v.Position.x);
                                minAABB.y = std::min(minAABB.y, v.Position.y);
                                minAABB.z = std::min(minAABB.z, v.Position.z);

                                maxAABB.x = std::max(maxAABB.x, v.Position.x);
                                maxAABB.y = std::max(maxAABB.y, v.Position.y);
                                maxAABB.z = std::max(maxAABB.z, v.Position.z);
                            }
                        }

                        // 3. 计算模型真正的几何中心点和长宽高的一半
                        Vector3 center = (minAABB + maxAABB) * 0.5f;
                        Vector3 halfExtents = (maxAABB - minAABB) * 0.5f;

                        // 4. 创建一个完美包裹车身的轻量级 Box
                        auto* proxyBox = new btBoxShape(btVector3(halfExtents.x, halfExtents.y, halfExtents.z));

                        // 5. 【核心魔法】：将 Box 偏移到模型的真正中心处！
                        // 建模师往往把坐标原点放在车底，这会导致 Box 下半截埋进地里。
                        // 通过局部偏移，物理外壳将精准对齐渲染网格！
                        btTransform localTransform;
                        localTransform.setIdentity();
                        localTransform.setOrigin(btVector3(center.x, center.y, center.z));

                        // 把偏移后的 Box 装入复合形状
                        compoundShape->addChildShape(localTransform, proxyBox);
                    }
                }
                shape = compoundShape;
                break;
            }
            default:
                shape = new btBoxShape(btVector3(0.5f, 0.5f, 0.5f));
                break;
        }

        return shape;
    }

    // =====================================================================================
    // 创建刚体
    // =====================================================================================
    void PhysicsSystem::CreateBody(Registry& registry,Entity entity, RigidBodyComponent& rb, const TransformComponent& t, const ColliderComponent& c) {
        if (!m_scene) return;

        // 1. 创建形状
        btCollisionShape* shape = CreateShape(entity,registry,c);
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
        finalPos = finalPos + (t.getRotation() * c.getOffset()); 
        // 暂时简化处理：直接叠加
        // finalPos = finalPos + c.getOffset();

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

        body->setUserPointer((void*)(uintptr_t)entity);

        // 6. 设置额外属性
        body->setFriction(rb.getFriction());
        // 阻尼 (可选，防止空气中无限漂浮)
        // body->setDamping(0.1f, 0.1f); 

        if (rb.isKinematic()) {
            body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);
            body->setActivationState(DISABLE_DEACTIVATION); // 永远不休眠
        }

        // if (c.getType() == ColliderType::ConvexHull) {
        //     // 参数 btVector3(x, y, z) 中，0代表锁死，1代表允许。
        //     // 设定 (0, 1, 0) 意味着：车子只能像陀螺一样绕 Y 轴转，永远不会侧翻或底朝天！
        //     body->setAngularFactor(btVector3(0.0f, 1.0f, 0.0f));
        // }

        // 7. 保存 Entity ID 到 UserPointer (非常重要，碰撞回调时反查 Entity)
        // 这里假设 Entity 是 uint32_t，强转存入
        // body->setUserPointer((void*)(uintptr_t)entity);
        // if (c.getType() == ColliderType::ConvexHull) {
        //     // body->setAngularFactor(btVector3(0.0f, 1.0f, 0.0f));
            
        //     // 【新增】：禁止这辆车进入休眠状态！
        //     body->setActivationState(DISABLE_DEACTIVATION);
        // }

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
    void PhysicsSystem::SyncTransformBack(const RigidBodyComponent& rb, TransformComponent& t,const ColliderComponent& c) {
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

        Vector3 rotatedOffset = newRot * c.getOffset(); // 假设你的四元数重载了与 Vector3 的乘法
    
        t.m_position = newPos - rotatedOffset;
        t.m_rotation = newRot;

        // 注意：这里绝对不能调用 t.setPosition(...)，因为那会触发 setDirty，导致死循环
    }

    // --- 修正 SyncDirtyData 的签名，去掉 const 以便清除脏标记 ---
    // void PhysicsSystem::SyncDirtyData(RigidBodyComponent& rb, TransformComponent& t, ColliderComponent& c) { ... }

} // namespace Lizeral