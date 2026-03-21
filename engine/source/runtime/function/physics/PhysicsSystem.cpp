#include "PhysicsSystem.h"

#include "runtime/function/physics/PhysicsScene.h"

#include "runtime/function/ecs/components/componentAll.h" 

namespace Lizeral {

    PhysicsSystem::PhysicsSystem() {}

    PhysicsSystem::~PhysicsSystem() {
        Shutdown();
    }

    void PhysicsSystem::Initialize(PhysicsScene* scene) {
        m_scene = scene;
        if (m_scene && m_scene->getWorld()) {
            btContactSolverInfo& info = m_scene->getWorld()->getSolverInfo();
            info.m_numIterations = 20; 
        }
    }

    void PhysicsSystem::Shutdown() {
        m_scene = nullptr; 
    }

    void PhysicsSystem::ResetPhysicsWorld(){
        if(m_scene){
            m_scene->Initialize();
        }
    }

    btVector3 PhysicsSystem::ToBtVector3(const Vector3& v) { return btVector3(v.x, v.y, v.z); }
    btQuaternion PhysicsSystem::ToBtQuaternion(const Quaternion& q) { return btQuaternion(q.x, q.y, q.z, q.w); }

    void PhysicsSystem::Tick(float deltaTime, Registry& registry) {
        if (!m_scene || !m_scene->getWorld()) return;

        btDiscreteDynamicsWorld* world = m_scene->getWorld();
        auto view = registry.view<RigidBodyComponent, TransformComponent, ColliderComponent>();

        for (auto entity : view) {
            auto& rb = view.get<RigidBodyComponent>(entity);
            auto& trans = view.get<TransformComponent>(entity);
            auto& col = view.get<ColliderComponent>(entity);

            if (!rb.getRuntimeBody()) {
                CreateBody(registry, entity, rb, trans, col);
                rb.clearDirty(0xFFFFFFFF); 
                trans.clearDirty(0xFFFFFFFF);
                col.clearDirty(0xFFFFFFFF);
                continue;
            }

            btRigidBody* body = static_cast<btRigidBody*>(rb.getRuntimeBody());
            if (body) {
                int flags = body->getCollisionFlags();
                if (col.getShowDebug()) {
                    body->setCollisionFlags(flags & ~btCollisionObject::CF_DISABLE_VISUALIZE_OBJECT);
                } else {
                    body->setCollisionFlags(flags | btCollisionObject::CF_DISABLE_VISUALIZE_OBJECT);
                }
            }

            if (rb.isDirty() || trans.isDirty() || col.isDirty()) {
                SyncDirtyData(rb, trans, col);
            }
        }

        world->stepSimulation(deltaTime, 10, 1.0f / 60.0f);

        for (auto entity : view) {
            auto& rb = view.get<RigidBodyComponent>(entity);
            if (rb.getRuntimeBody() && rb.getMass() > 0.0f && !rb.isKinematic()) {
                auto& trans = view.get<TransformComponent>(entity);
                auto& col = view.get<ColliderComponent>(entity);
                SyncTransformBack(rb, trans, col);
            }
        }
    }

    void PhysicsSystem::UpdateDebugLines(Registry& registry) {
        m_debugLines.clear(); 

        auto view = registry.view<TransformComponent, ColliderComponent>();
        for (auto entity : view) {
            auto& col = view.get<ColliderComponent>(entity);

            if (!col.getShowDebug()) continue;

            auto& trans = view.get<TransformComponent>(entity);
            
            if (col.getType() == ColliderType::Box) {
                Vector3 center = trans.getPosition() + (trans.getRotation() * col.getOffset());
                Vector3 extents = col.getSize() * 0.5f * trans.getScale(); 
                Quaternion rot = trans.getRotation();
                Vector3 color(0.0f, 1.0f, 0.0f);

                Vector3 v[8];
                int idx = 0;
                for (int x = -1; x <= 1; x += 2) {
                    for (int y = -1; y <= 1; y += 2) {
                        for (int z = -1; z <= 1; z += 2) {
                            Vector3 local(x * extents.x, y * extents.y, z * extents.z);
                            v[idx++] = center + (rot * local);
                        }
                    }
                }

                auto AddLine = [&](const Vector3& p1, const Vector3& p2) {
                    m_debugLines.push_back({p1, color});
                    m_debugLines.push_back({p2, color});
                };
                AddLine(v[0], v[1]); AddLine(v[1], v[3]); AddLine(v[3], v[2]); AddLine(v[2], v[0]); // 底面
                AddLine(v[4], v[5]); AddLine(v[5], v[7]); AddLine(v[7], v[6]); AddLine(v[6], v[4]); // 顶面
                AddLine(v[0], v[4]); AddLine(v[1], v[5]); AddLine(v[2], v[6]); AddLine(v[3], v[7]); // 侧边
            }
        }
    }

    void PhysicsSystem::DebugDrawWorld(Registry& registry) {
        if (!m_scene || !m_scene->getWorld() ) return;

        auto view = registry.view<RigidBodyComponent, ColliderComponent>();
        for (auto entity : view) {
            auto& rb = view.get<RigidBodyComponent>(entity);
            auto& col = view.get<ColliderComponent>(entity);
            
            btRigidBody* body = static_cast<btRigidBody*>(rb.getRuntimeBody());
            if (body) {
                int flags = body->getCollisionFlags();
                if (col.getShowDebug()) {
                    body->setCollisionFlags(flags & ~btCollisionObject::CF_DISABLE_VISUALIZE_OBJECT);
                } else {
                    body->setCollisionFlags(flags | btCollisionObject::CF_DISABLE_VISUALIZE_OBJECT);
                }
            }
        }

        m_scene->getWorld()->debugDrawWorld();
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

            const btCollisionObject* obj = callback.m_collisionObject;
            outHit.entity = static_cast<uint32_t>((uintptr_t)obj->getUserPointer()); 
            
            return true;
        }

        outHit.hasHit = false;
        return false;
    }

    btCollisionShape* PhysicsSystem::CreateShape(const ColliderComponent& col) {
        btCollisionShape* shape = nullptr;
        switch (col.getType()) {
            case ColliderType::Box: 
                shape = new btBoxShape(ToBtVector3(col.getSize() * 0.5f)); 
                shape->setMargin(0.005f);
                break;
            case ColliderType::Sphere: 
                shape = new btSphereShape(col.getRadius()); 
                break;
            case ColliderType::Capsule: 
                shape = new btCapsuleShape(col.getRadius(), col.getHeight()); 
                break;
            default: 
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
                shape->setMargin(0.005f);
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
            default:
                shape = new btBoxShape(btVector3(0.5f, 0.5f, 0.5f));
                break;
        }

        return shape;
    }

    void PhysicsSystem::CreateBody(Registry& registry, Entity entity, RigidBodyComponent& rb, const TransformComponent& t, const ColliderComponent& c) {
        if (!m_scene) return;

        btCollisionShape* shape = CreateShape(entity, registry, c);
        m_scene->trackShape(shape);

        btTransform startTrans;
        startTrans.setIdentity();
        Vector3 finalPos = t.getPosition() + (t.getRotation() * c.getOffset()); 
        startTrans.setOrigin(ToBtVector3(finalPos));
        startTrans.setRotation(ToBtQuaternion(t.getRotation()));

        shape->setLocalScaling(ToBtVector3(t.getScale()));

        btVector3 localInertia(0, 0, 0);
        float mass = rb.getMass();
        if (mass > 0.0f) shape->calculateLocalInertia(mass, localInertia);

        btDefaultMotionState* myMotionState = new btDefaultMotionState(startTrans);
        btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, shape, localInertia);
        btRigidBody* body = new btRigidBody(rbInfo);

        body->setUserPointer((void*)(uintptr_t)entity);
        body->setFriction(rb.getFriction());
        body->setRestitution(rb.getRestitution());

        if (rb.isKinematic()) {
            body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);
            body->setActivationState(DISABLE_DEACTIVATION);
        }

        if (c.getShowDebug()) {
            body->setCollisionFlags(body->getCollisionFlags() & ~btCollisionObject::CF_DISABLE_VISUALIZE_OBJECT);
        } else {
            body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_DISABLE_VISUALIZE_OBJECT);
        }

        m_scene->getWorld()->addRigidBody(body);
        rb.setRuntimeBody(body);
    }

    void PhysicsSystem::SyncDirtyData(RigidBodyComponent& rb, TransformComponent& t, ColliderComponent& c) {
        btRigidBody* body = static_cast<btRigidBody*>(rb.getRuntimeBody());
        if (!body) return;

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
            body->activate(true);
            rb.clearDirty(PHYS_DIRTY_RESTITUTION);
        }
        
        // Transform
        if (t.isDirty(TRANS_DIRTY_POSITION) || t.isDirty(TRANS_DIRTY_ROTATION) || c.isDirty()) {
            btTransform tr;
            Vector3 finalPos = t.getPosition() + (t.getRotation() * c.getOffset());
            tr.setOrigin(ToBtVector3(finalPos));
            tr.setRotation(ToBtQuaternion(t.getRotation()));
            body->setWorldTransform(tr);
            if (body->getMotionState()) body->getMotionState()->setWorldTransform(tr);
            body->activate(true);
        }

        // Scale 
        if (t.isDirty(TRANS_DIRTY_SCALE)) {
            if (body->getCollisionShape()) {
                body->getCollisionShape()->setLocalScaling(ToBtVector3(t.getScale()));
                if (rb.getMass() > 0.0f) {
                    btVector3 inertia;
                    body->getCollisionShape()->calculateLocalInertia(rb.getMass(), inertia);
                    body->setMassProps(rb.getMass(), inertia);
                    body->updateInertiaTensor();
                }
                body->activate(true);
            }
        }

        if (c.isDirty()) {
            btCollisionShape* oldShape = body->getCollisionShape();
            btCollisionShape* newShape = CreateShape(c); 
            newShape->setLocalScaling(oldShape->getLocalScaling());
            
            body->setCollisionShape(newShape);
            m_scene->trackShape(newShape); 

            if (rb.getMass() > 0.0f) {
                btVector3 inertia(0, 0, 0);
                newShape->calculateLocalInertia(rb.getMass(), inertia);
                body->setMassProps(rb.getMass(), inertia);
                body->updateInertiaTensor();
            }
            body->activate(true);

        }

        rb.clearDirty(0xFFFFFFFF);
        t.clearDirty(0xFFFFFFFF);
        c.clearDirty(0xFFFFFFFF);
    }

    void PhysicsSystem::SyncTransformBack(const RigidBodyComponent& rb, TransformComponent& t, const ColliderComponent& c) {
        btRigidBody* body = static_cast<btRigidBody*>(rb.getRuntimeBody());
        btTransform trans;
        if (body->getMotionState()) {
            body->getMotionState()->getWorldTransform(trans);
        } else {
            trans = body->getWorldTransform();
        }

        btVector3 pos = trans.getOrigin();
        btQuaternion rot = trans.getRotation();

        Quaternion newRot(rot.w(), rot.x(), rot.y(), rot.z());
        Vector3 newPos(pos.x(), pos.y(), pos.z());

        Vector3 rotatedOffset = newRot * c.getOffset(); 

        t.m_position = newPos - rotatedOffset;
        t.m_rotation = newRot;
    }

} // namespace Lizeral