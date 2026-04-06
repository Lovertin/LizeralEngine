#include "PhysicsSystem.h"

#include <chrono>

#include "runtime/function/physics/PhysicsScene.h"
#include "runtime/function/ecs/components/componentAll.h"
#include "runtime/function/ecs/hybrid/hybrid_registry.h"
#include "runtime/function/ecs/hybrid/hybrid_default_traits.h"

namespace {

    using PhysicsClock = std::chrono::high_resolution_clock;

    double ToMicroseconds(const PhysicsClock::duration& duration) {
        return std::chrono::duration<double, std::micro>(duration).count();
    }

} // namespace

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
        m_trackedBodies.clear();
        m_scene = nullptr;
    }

    void PhysicsSystem::ResetPhysicsWorld() {
        if (m_scene) {
            m_scene->Initialize();
        }
        m_trackedBodies.clear();
    }

    void PhysicsSystem::ResetPhysicsWorld(Registry& registry) {
        ResetPhysicsWorldInternal(registry);
    }

    void PhysicsSystem::ResetPhysicsWorld(HybridRegistry& registry) {
        ResetPhysicsWorldInternal(registry);
    }

    btVector3 PhysicsSystem::ToBtVector3(const Vector3& v) { return btVector3(v.x, v.y, v.z); }
    btQuaternion PhysicsSystem::ToBtQuaternion(const Quaternion& q) { return btQuaternion(q.x, q.y, q.z, q.w); }

    template<typename RegistryT>
    PhysicsAccessProfile PhysicsSystem::ProfileComponentAccessInternal(RegistryT& registry) {
        PhysicsAccessProfile profile{};

        const auto totalStart = PhysicsClock::now();
        const auto viewStart = PhysicsClock::now();
        auto view = registry.template view<RigidBodyComponent, TransformComponent, ColliderComponent>();
        const auto iterateStart = PhysicsClock::now();

        volatile double sink = 0.0;
        for (auto entity : view) {
            auto& rb = view.template get<RigidBodyComponent>(entity);
            auto& trans = view.template get<TransformComponent>(entity);
            auto& col = view.template get<ColliderComponent>(entity);

            sink += rb.getMass();
            sink += trans.getPosition().x + trans.getScale().x;
            sink += static_cast<double>(static_cast<int>(col.getType()));
            ++profile.entity_count;
        }

        profile.view_build_us = ToMicroseconds(iterateStart - viewStart);
        profile.iteration_us = ToMicroseconds(PhysicsClock::now() - iterateStart);
        profile.total_us = ToMicroseconds(PhysicsClock::now() - totalStart);
        (void)sink;
        return profile;
    }

    PhysicsAccessProfile PhysicsSystem::ProfileChunkAccessInternal(HybridRegistry& registry) {
        PhysicsAccessProfile profile{};

        const auto totalStart = PhysicsClock::now();
        const auto viewStart = PhysicsClock::now();
        auto chunks = registry.chunk_query<TransformComponent, RigidBodyComponent, ColliderComponent>();
        const auto iterateStart = PhysicsClock::now();

        volatile double sink = 0.0;
        for (auto& chunk : chunks) {
            for (size_t i = 0; i < chunk.size(); ++i) {
                auto& trans = chunk.template get<TransformComponent>(i);
                auto& rb = chunk.template get<RigidBodyComponent>(i);
                auto& col = chunk.template get<ColliderComponent>(i);

                sink += rb.getMass();
                sink += trans.getPosition().x + trans.getScale().x;
                sink += static_cast<double>(static_cast<int>(col.getType()));
                ++profile.entity_count;
            }
        }

        profile.view_build_us = ToMicroseconds(iterateStart - viewStart);
        profile.iteration_us = ToMicroseconds(PhysicsClock::now() - iterateStart);
        profile.total_us = ToMicroseconds(PhysicsClock::now() - totalStart);
        (void)sink;
        return profile;
    }

    template<typename RegistryT>
    PhysicsTickProfile PhysicsSystem::TickInternal(float deltaTime, RegistryT& registry, PhysicsTickProfile* profileSink) {
        PhysicsTickProfile profile{};
        if (!m_scene || !m_scene->getWorld()) {
            if (profileSink) {
                *profileSink = profile;
            }
            return profile;
        }

        btDiscreteDynamicsWorld* world = m_scene->getWorld();
        const bool captureProfile = profileSink != nullptr;

        PhysicsClock::time_point totalStart{};
        PhysicsClock::time_point preStepStart{};
        PhysicsClock::time_point simulationStart{};
        PhysicsClock::time_point writeBackStart{};

        if (captureProfile) {
            totalStart = PhysicsClock::now();
            preStepStart = totalStart;
        }

        CleanupStaleBodiesInternal(registry);

        auto view = registry.template view<RigidBodyComponent, TransformComponent, ColliderComponent>();

        for (auto entity : view) {
            auto& rb = view.template get<RigidBodyComponent>(entity);
            auto& trans = view.template get<TransformComponent>(entity);
            auto& col = view.template get<ColliderComponent>(entity);

            ++profile.entity_count;

            if (!rb.getRuntimeBody()) {
                CreateBody(entity, rb, trans, col);
                rb.clearDirty(0xFFFFFFFF);
                trans.clearDirty(0xFFFFFFFF);
                col.clearDirty(0xFFFFFFFF);
                ++profile.created_body_count;
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
                ++profile.dirty_sync_count;
            }
        }

        if (captureProfile) {
            simulationStart = PhysicsClock::now();
            profile.pre_step_us = ToMicroseconds(simulationStart - preStepStart);
        }

        world->stepSimulation(deltaTime, 10, 1.0f / 60.0f);

        if (captureProfile) {
            writeBackStart = PhysicsClock::now();
            profile.simulation_us = ToMicroseconds(writeBackStart - simulationStart);
        }

        for (auto entity : view) {
            auto& rb = view.template get<RigidBodyComponent>(entity);
            if (rb.getRuntimeBody() && rb.getMass() > 0.0f && !rb.isKinematic()) {
                auto& trans = view.template get<TransformComponent>(entity);
                auto& col = view.template get<ColliderComponent>(entity);
                SyncTransformBack(rb, trans, col);
                ++profile.transform_write_back_count;
            }
        }

        if (captureProfile) {
            const auto totalEnd = PhysicsClock::now();
            profile.write_back_us = ToMicroseconds(totalEnd - writeBackStart);
            profile.total_us = ToMicroseconds(totalEnd - totalStart);
            *profileSink = profile;
        }

        return profile;
    }

    PhysicsTickProfile PhysicsSystem::TickChunkedInternal(float deltaTime, HybridRegistry& registry, PhysicsTickProfile* profileSink) {
        PhysicsTickProfile profile{};
        if (!m_scene || !m_scene->getWorld()) {
            if (profileSink) {
                *profileSink = profile;
            }
            return profile;
        }

        btDiscreteDynamicsWorld* world = m_scene->getWorld();
        const bool captureProfile = profileSink != nullptr;

        PhysicsClock::time_point totalStart{};
        PhysicsClock::time_point preStepStart{};
        PhysicsClock::time_point simulationStart{};
        PhysicsClock::time_point writeBackStart{};

        if (captureProfile) {
            totalStart = PhysicsClock::now();
            preStepStart = totalStart;
        }

        CleanupStaleBodiesInternal(registry);

        auto chunks = registry.chunk_query<TransformComponent, RigidBodyComponent, ColliderComponent>();

        for (auto& chunk : chunks) {
            profile.entity_count += chunk.size();

            for (size_t i = 0; i < chunk.size(); ++i) {
                const Entity entity = chunk.entities()[i];
                auto& trans = chunk.template get<TransformComponent>(i);
                auto& rb = chunk.template get<RigidBodyComponent>(i);
                auto& col = chunk.template get<ColliderComponent>(i);

                if (!rb.getRuntimeBody()) {
                    CreateBody(entity, rb, trans, col);
                    rb.clearDirty(0xFFFFFFFF);
                    trans.clearDirty(0xFFFFFFFF);
                    col.clearDirty(0xFFFFFFFF);
                    ++profile.created_body_count;
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
                    ++profile.dirty_sync_count;
                }
            }
        }

        if (captureProfile) {
            simulationStart = PhysicsClock::now();
            profile.pre_step_us = ToMicroseconds(simulationStart - preStepStart);
        }

        world->stepSimulation(deltaTime, 10, 1.0f / 60.0f);

        if (captureProfile) {
            writeBackStart = PhysicsClock::now();
            profile.simulation_us = ToMicroseconds(writeBackStart - simulationStart);
        }

        for (auto& chunk : chunks) {
            for (size_t i = 0; i < chunk.size(); ++i) {
                auto& trans = chunk.template get<TransformComponent>(i);
                auto& rb = chunk.template get<RigidBodyComponent>(i);
                auto& col = chunk.template get<ColliderComponent>(i);

                if (rb.getRuntimeBody() && rb.getMass() > 0.0f && !rb.isKinematic()) {
                    SyncTransformBack(rb, trans, col);
                    ++profile.transform_write_back_count;
                }
            }
        }

        if (captureProfile) {
            const auto totalEnd = PhysicsClock::now();
            profile.write_back_us = ToMicroseconds(totalEnd - writeBackStart);
            profile.total_us = ToMicroseconds(totalEnd - totalStart);
            *profileSink = profile;
        }

        return profile;
    }

    void PhysicsSystem::Tick(float deltaTime, Registry& registry) {
        TickInternal(deltaTime, registry, nullptr);
    }

    void PhysicsSystem::Tick(float deltaTime, HybridRegistry& registry) {
        if (m_hybridExecutionMode == PhysicsHybridExecutionMode::Chunked) {
            TickChunkedInternal(deltaTime, registry, nullptr);
            return;
        }

        TickInternal(deltaTime, registry, nullptr);
    }

    void PhysicsSystem::TickChunked(float deltaTime, HybridRegistry& registry) {
        TickChunkedInternal(deltaTime, registry, nullptr);
    }

    PhysicsTickProfile PhysicsSystem::TickProfiled(float deltaTime, Registry& registry) {
        return TickInternal(deltaTime, registry, &m_lastSparseTickProfile);
    }

    PhysicsTickProfile PhysicsSystem::TickProfiled(float deltaTime, HybridRegistry& registry) {
        if (m_hybridExecutionMode == PhysicsHybridExecutionMode::Chunked) {
            return TickChunkedInternal(deltaTime, registry, &m_lastHybridChunkTickProfile);
        }

        return TickInternal(deltaTime, registry, &m_lastHybridTickProfile);
    }

    PhysicsTickProfile PhysicsSystem::TickProfiledChunked(float deltaTime, HybridRegistry& registry) {
        return TickChunkedInternal(deltaTime, registry, &m_lastHybridChunkTickProfile);
    }

    PhysicsAccessProfile PhysicsSystem::ProfileComponentAccess(Registry& registry) {
        return ProfileComponentAccessInternal(registry);
    }

    PhysicsAccessProfile PhysicsSystem::ProfileComponentAccess(HybridRegistry& registry) {
        return ProfileComponentAccessInternal(registry);
    }

    PhysicsAccessProfile PhysicsSystem::ProfileChunkAccess(HybridRegistry& registry) {
        return ProfileChunkAccessInternal(registry);
    }

    PhysicsAccessComparison PhysicsSystem::CompareComponentAccess(Registry& sparseRegistry, HybridRegistry& hybridRegistry) {
        PhysicsAccessComparison comparison{};
        comparison.sparse_set = ProfileComponentAccess(sparseRegistry);
        comparison.hybrid_archetype = ProfileComponentAccess(hybridRegistry);
        return comparison;
    }

    template<typename RegistryT>
    void PhysicsSystem::UpdateDebugLinesInternal(RegistryT& registry) {
        m_debugLines.clear();

        auto view = registry.template view<TransformComponent, ColliderComponent>();
        for (auto entity : view) {
            auto& col = view.template get<ColliderComponent>(entity);

            if (!col.getShowDebug()) continue;

            auto& trans = view.template get<TransformComponent>(entity);

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

                auto addLine = [&](const Vector3& p1, const Vector3& p2) {
                    m_debugLines.push_back({p1, color});
                    m_debugLines.push_back({p2, color});
                };

                addLine(v[0], v[1]); addLine(v[1], v[3]); addLine(v[3], v[2]); addLine(v[2], v[0]);
                addLine(v[4], v[5]); addLine(v[5], v[7]); addLine(v[7], v[6]); addLine(v[6], v[4]);
                addLine(v[0], v[4]); addLine(v[1], v[5]); addLine(v[2], v[6]); addLine(v[3], v[7]);
            }
        }
    }

    void PhysicsSystem::UpdateDebugLines(Registry& registry) {
        UpdateDebugLinesInternal(registry);
    }

    void PhysicsSystem::UpdateDebugLines(HybridRegistry& registry) {
        UpdateDebugLinesInternal(registry);
    }

    template<typename RegistryT>
    void PhysicsSystem::DebugDrawWorldInternal(RegistryT& registry) {
        if (!m_scene || !m_scene->getWorld()) return;

        auto view = registry.template view<RigidBodyComponent, ColliderComponent>();
        for (auto entity : view) {
            auto& rb = view.template get<RigidBodyComponent>(entity);
            auto& col = view.template get<ColliderComponent>(entity);

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

    void PhysicsSystem::DebugDrawWorld(Registry& registry) {
        DebugDrawWorldInternal(registry);
    }

    void PhysicsSystem::DebugDrawWorld(HybridRegistry& registry) {
        DebugDrawWorldInternal(registry);
    }

    bool PhysicsSystem::Raycast(const Ray& ray, RaycastHit& outHit, float maxDistance) {
        if (!m_scene || !m_scene->getWorld()) return false;

        Vector3 start = ray.origin;
        Vector3 end = ray.GetPoint(maxDistance);

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

    void PhysicsSystem::Clear() {
        m_debugLines.clear();
        m_lastSparseTickProfile = {};
        m_lastHybridTickProfile = {};
        m_lastHybridChunkTickProfile = {};
    }

    template<typename RegistryT>
    void PhysicsSystem::CleanupStaleBodiesInternal(RegistryT& registry) {
        if (!m_scene || !m_scene->getWorld()) {
            m_trackedBodies.clear();
            return;
        }

        for (auto it = m_trackedBodies.begin(); it != m_trackedBodies.end();) {
            const Entity entity = it->first;
            btRigidBody* trackedBody = it->second;

            const bool hasPhysicsComponents =
                registry.template has<RigidBodyComponent>(entity) &&
                registry.template has<TransformComponent>(entity) &&
                registry.template has<ColliderComponent>(entity);

            bool shouldDestroy = !hasPhysicsComponents;
            RigidBodyComponent* rigidBody = nullptr;

            if (hasPhysicsComponents) {
                rigidBody = &registry.template get<RigidBodyComponent>(entity);
                shouldDestroy = static_cast<btRigidBody*>(rigidBody->getRuntimeBody()) != trackedBody;
            }

            if (!shouldDestroy) {
                ++it;
                continue;
            }

            if (rigidBody != nullptr && rigidBody->getRuntimeBody() == trackedBody) {
                rigidBody->setRuntimeBody(nullptr);
            }

            DestroyBody(trackedBody);
            it = m_trackedBodies.erase(it);
        }
    }

    template<typename RegistryT>
    void PhysicsSystem::ResetPhysicsWorldInternal(RegistryT& registry) {
        for (auto& pair : m_trackedBodies) {
            const Entity entity = pair.first;
            if (registry.template has<RigidBodyComponent>(entity)) {
                registry.template get<RigidBodyComponent>(entity).setRuntimeBody(nullptr);
            }
        }

        m_trackedBodies.clear();
        ResetPhysicsWorld();
    }

    void PhysicsSystem::DestroyBody(btRigidBody* body) {
        if (!body) {
            return;
        }

        if (m_scene && m_scene->getWorld()) {
            m_scene->getWorld()->removeRigidBody(body);
        }

        if (body->getMotionState()) {
            delete body->getMotionState();
        }

        delete body;
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

    void PhysicsSystem::CreateBody(Entity entity, RigidBodyComponent& rb, const TransformComponent& t, const ColliderComponent& c) {
        if (!m_scene) return;

        btCollisionShape* shape = CreateShape(c);
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
        m_trackedBodies[entity] = body;
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

        if (t.isDirty(TRANS_DIRTY_POSITION) || t.isDirty(TRANS_DIRTY_ROTATION) || c.isDirty()) {
            btTransform tr;
            Vector3 finalPos = t.getPosition() + (t.getRotation() * c.getOffset());
            tr.setOrigin(ToBtVector3(finalPos));
            tr.setRotation(ToBtQuaternion(t.getRotation()));
            body->setWorldTransform(tr);
            if (body->getMotionState()) body->getMotionState()->setWorldTransform(tr);
            body->activate(true);
        }

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
